import os, json, datetime
import rethinkdb as r
from Bio import SeqIO
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='tdb', help="database to download from")
parser.add_argument('-v', '--virus', default='h3n2', help="virus table to interact with")
parser.add_argument('--subtype', nargs='+', type=str, default=None, help="subtype to be include in download")
parser.add_argument('--host', nargs='+', type=str, default=None, help="hosts to be include in download")
parser.add_argument('--path', default='data', help="path to dump output files to")
parser.add_argument('--ftype', default='text', help="output file format, default \"text\", other is \"json\"")
parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class tdb_download(object):
    def __init__(self, database, virus, rethink_host=None, auth_key=None, **kwargs):
        '''
        parser for virus, fasta fields, output file names, output file format path, interval
        '''
        self.virus = virus.lower()
        self.database = database.lower()
        if self.database not in ['tdb', 'test_tdb']:
            raise Exception("Cant download to this database: " + self.database)
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
        self.rethink_io.connect_rethink(self.database, self.auth_key, self.rethink_host)
        self.rethink_io.check_table_exists(self.database, self.virus)
        self.measurements = []

    def count_documents(self):
        '''
        return integer count of number of documents in table
        '''
        return r.db(self.database).table(self.virus).count().run()

    def download(self, output=True, **kwargs):
        '''
        download documents from table
        '''
        print("Downloading all titer measurements from the table: " + self.virus)
        self.measurements = list(r.db(self.database).table(self.virus).run())
        self.measurements = self.subsetting(self.measurements, **kwargs)
        if output:
            self.output(**kwargs)

    def subsetting(self, cursor, subtype=None, host=None, **kwargs):
        '''
        filter through documents in vdb to return subsets of sequence
        '''
        result = cursor
        print("Documents in table before subsetting: " + str(len(result)))
        if subtype is not None:
            subtype = [s.lower() for s in subtype]
            result = filter(lambda doc: doc['subtype'] in subtype, result)
            print('Removed documents that were not of the specified subtype (' + ','.join(subtype) + '), remaining documents: ' + str(len(result)))
        if host is not None:
            host = [h.title() for h in host]
            result = filter(lambda doc: doc['host'] in host, result)
            print('Removed documents that were not of the specified host (' + ','.join(host) + '), remaining documents: ' + str(len(result)))
        print("Documents in table after subsetting: " + str(len(result)))
        return result

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

    def write_text(self, measurements, fname, text_fields=['virus', 'serum', 'ferret_id', 'source', 'titer']):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for meas in measurements:
                for field in text_fields:
                    if field in meas and meas[field] is not None:
                        handle.write(str(meas[field]) + '\t')
                    else:
                        handle.write('?' + '\t')
                handle.write('\n')
            handle.close()
            print("Wrote to " + fname)

    def output(self, path, fstem, ftype, **kwargs):
        fname = path + '/' + fstem + '.' + ftype
        if ftype == 'json':
            self.write_json(self.measurements,fname)
        elif ftype == 'text':
            self.write_text(self.measurements, fname)
        else:
            raise Exception("Can't output to that file type, only json or text allowed")

if __name__=="__main__":
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = tdb_download(**args.__dict__)
    connTDB.download(**args.__dict__)
