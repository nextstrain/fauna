import os, json, datetime
import rethinkdb as r
from Bio import SeqIO
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to download from")
parser.add_argument('-v', '--virus', default='zika', help="virus table to interact with")
parser.add_argument('--path', default='data', help="path to dump output files to")
parser.add_argument('--ftype', default='fasta', help="output file format, default \"fasta\", other is \"json\"")
parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--public_only', default=False, action="store_true", help="include to subset public sequences")
parser.add_argument('--countries', nargs='+', type=str, default=None, help="Countries(in CamelCase Format) to be include in download")

class vdb_download(object):

    def __init__(self, database, virus, rethink_host=None, auth_key=None, **kwargs):
        '''
        parser for virus, fasta fields, output file names, output file format path, interval
        '''
        self.virus = virus.lower()
        self.database = database.lower()
        if self.database not in ['vdb', 'test_vdb']:
            raise Exception("Cant download from this database: " + self.database)
        if rethink_host is None:
            try:
                self.rethink_host = os.environ['RETHINK_HOST']
            except:
                raise Exception("Missing rethink host")
        else:
            self.rethink_host = rethink_host
        if auth_key is None:
            try:
                self.auth_key = os.environ['RETHINK_AUTH_KEY']
            except:
                raise Exception("Missing rethink auth_key")
        else:
            self.auth_key = auth_key
        self.connect_rethink()
        self.viruses = []

    def connect_rethink(self):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        try:
            r.connect(host=self.rethink_host, port=28015, db=self.database, auth_key=self.auth_key).repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

        existing_tables = r.db(self.database).table_list().run()
        if self.virus not in existing_tables:
            raise Exception("No table exists yet for " + self.virus)

    def count_documents(self):
        '''
        return integer count of number of documents in table
        '''
        return r.db(self.database).table(self.virus).count().run()

    def download(self, output=True, **kwargs):
        '''
        download documents from table
        '''
        print("Downloading all viruses from the table: " + self.virus)
        self.viruses = list(r.db(self.database).table(self.virus).run())
        for doc in self.viruses:
            self.pick_best_sequence(doc)
        self.viruses = self.subsetting(self.viruses, **kwargs)
        if output:
            self.output(**kwargs)

    def subsetting(self, cursor, public_only=False, countries=None, **kwargs):
        '''
        filter through documents in vdb to return subsets of sequence
        '''
        result = cursor
        print("Documents in table before subsetting: " + str(len(result)))
        if public_only:
            result = filter(lambda doc: doc['public'], result)
            print('Removed documents that were not public, remaining documents: ' + str(len(result)))
        if countries is not None:
            result = filter(lambda doc: doc['country'] in countries, result)
            print('Removed documents that were not in countries specified ('
                + ','.join(countries) + '), remaining documents: ' + str(len(result)))
        print("Documents in table after subsetting: " + str(len(result)))
        return result

    def pick_best_sequence(self, document):
        '''
        find the best sequence in the given document. Currently by longest sequence.
        Resulting document is with flatter dictionary structure
        '''
        longest_sequence_pos = np.argmax([len(seq_info['sequence']) for seq_info in document['sequences']])
        best_sequence_info = document['sequences'][longest_sequence_pos]
        best_citation_info = document['citations'][longest_sequence_pos]

        # create flatter structure for virus info
        for atr in best_sequence_info.keys():
            document[atr] = best_sequence_info[atr]
        for atr in best_citation_info.keys():
            document[atr] = best_citation_info[atr]
        del document['sequences']
        del document['citations']

    def write_json(self, data, fname, indent=1):
        '''
        writes as list of viruses (dictionaries)
        '''
        try:
            handle = open(fname, 'w')
        except:
            print("Couldn't open output file")
            print(fname)
            raise FileNotFoundError
        else:
            json.dump(data, handle, indent=indent)
            handle.close()
            print("Wrote to " + fname)

    def write_fasta(self, viruses, fname, sep='|', fasta_fields=['strain', 'virus', 'accession', 'date', 'region',
                                                                 'country', 'division', 'location', 'source', 'locus',
                                                                 'authors', 'subtype']):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for v in viruses:
                fields = [str(v[field]) if (field in v and v[field] is not None) else '?'
                          for field in fasta_fields]
                handle.write(">"+sep.join(fields)+'\n')
                handle.write(v['sequence'] + "\n")
            handle.close()
            print("Wrote to " + fname)

    def output(self, path, fstem, ftype):
        fname = path + '/' + fstem + '.' + ftype
        if ftype == 'json':
            self.write_json(self.viruses,fname)
        elif ftype == 'fasta':
            self.write_fasta(self.viruses, fname)
        else:
            raise Exception("Can't output to that file type, only json or fasta allowed")

if __name__=="__main__":
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = vdb_download(**args.__dict__)
    connVDB.download(path=args.path, fstem=args.fstem, ftype=args.ftype)