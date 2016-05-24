import os, json, datetime
import rethinkdb as r
from Bio import SeqIO
import numpy as np
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to download from")
parser.add_argument('-v', '--virus', default='zika', help="virus table to interact with")
parser.add_argument('--path', default='data', help="path to dump output files to")
parser.add_argument('--ftype', default='fasta', help="output file format, default \"fasta\", other is \"json\"")
parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
parser.add_argument('--fasta_fields', default=['strain', 'virus', 'accession', 'date', 'region', 'authors'], help="fasta fields for output fasta")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--public_only', default=False, action="store_true", help="include to subset public sequences")
parser.add_argument('--select', nargs='+', type=str, default=None, help="Select specific fields ie \'--select field1:value1 field2:value1,value2\'")

class vdb_download(object):

    def __init__(self, database, virus, rethink_host=None, auth_key=None, **kwargs):
        '''
        parser for virus, fasta fields, output file names, output file format path, interval
        '''
        self.virus_specific_fasta_fields = []
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
        if self.rethink_host == "localhost":
            self.auth_key = None
        elif auth_key is not None:
            self.auth_key = auth_key
        else:
            try:
                self.auth_key = os.environ['RETHINK_AUTH_KEY']
            except:
                raise Exception("Missing rethink auth_key")
        self.rethink_io = rethink_io()
        self.rethink_io.connect_rethink(self.database, self.virus, self.auth_key, self.rethink_host)
        self.viruses = []

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

    def subsetting(self, cursor, public_only=False, select=None, **kwargs):
        '''
        filter through documents in vdb to return subsets of sequence
        '''
        # determine fields in all documents that can be filtered by
        list_of_fields = [set(doc.keys()) for doc in cursor]
        fields = set.intersection(*list_of_fields)

        print("Documents in table before subsetting: " + str(len(cursor)))
        if public_only:
            cursor = filter(lambda doc: doc['public'], cursor)
            print('Removed documents that were not public, remaining documents: ' + str(len(cursor)))
        if select is not None:
            selections = self.parse_grouping_argument(select)
            for sel in selections:
                if sel[0] in fields:
                    cursor = filter(lambda doc: doc[sel[0]] in sel[1], cursor)
                    print('Removed documents that were not in values specified (' + ','.join(sel[1]) + ') for field \'' + sel[0] + '\', remaining documents: ' + str(len(cursor)))
                else:
                    print(sel[0] + " is not in all documents, can't subset by that field")
        print("Documents in table after subsetting: " + str(len(cursor)))
        return cursor

    def parse_grouping_argument(self, grouping):
        '''
        parse the select parameter to determine which group name to filter and for what values
        :return: (grouping name, group values))
        '''
        selections = []
        for group in grouping:
            result = group.split(':')
            selections.append((result[0].lower(), result[1].lower().split(',')))
        return selections

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
                                                                 'authors']):
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

    def output(self, path, fstem, ftype, fasta_fields, **kwargs):
        fname = path + '/' + fstem + '.' + ftype
        if ftype == 'json':
            self.write_json(self.viruses,fname)
        elif ftype == 'fasta':
            self.write_fasta(self.viruses, fname, fasta_fields=fasta_fields+self.virus_specific_fasta_fields)
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
    connVDB.download(**args.__dict__)