import os, json, datetime
import rethinkdb as r
from Bio import SeqIO

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to download from")
parser.add_argument('-v', '--virus', default='Zika', help="virus table to interact with")
parser.add_argument('--path', default='tdb/data/', help="path to dump output files to")
parser.add_argument('--ftype', default='txt', help="output file format, default \"txt\", other is \"json\"")
parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
parser.add_argument('--host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class vdb_download(object):
    def __init__(self, **kwargs):
        '''
        parser for virus, fasta fields, output file names, output file format path, interval
        '''
        self.kwargs = kwargs
        if 'host' in self.kwargs:
            self.host = self.kwargs['host']
        if 'RETHINK_HOST' in os.environ and self.host is None:
            self.host = os.environ['RETHINK_HOST']
        if self.host is None:
            raise Exception("Missing rethink host")
        if 'auth_key' in self.kwargs:
            self.auth_key = self.kwargs['auth_key']
        if 'RETHINK_AUTH_KEY' in os.environ and self.auth_key is None:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing rethink auth_key")
        if 'path' in self.kwargs:
            self.path = self.kwargs['path']
        if not os.path.isdir(self.path):
            os.makedirs(self.path)

        if 'database' in self.kwargs:
            self.database = self.kwargs['database']
        if 'virus' in self.kwargs:
            self.virus = self.kwargs['virus'].lower()
        self.measurements = []

        if 'ftype' in self.kwargs:
            self.ftype = self.kwargs['ftype']
        self.current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
        if 'fstem' in self.kwargs:
            self.fstem = self.kwargs['fstem']
        if self.fstem is None:
            self.fstem = self.virus + '_' + self.current_date
        self.fname = self.fstem + '.' + self.ftype

        self.txt_fields = ['virus', 'serum', 'ferret_id', 'source', 'titer']
        self.connect_rethink()

    def connect_rethink(self):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        try:
            r.connect(host=self.host, port=28015, db=self.database, auth_key=self.auth_key).repl()
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

    def download(self):
        '''
        download documents from table
        '''
        print("Downloading all titer measurements from the table: " + self.virus)
        self.measurements = list(r.db(self.database).table(self.virus).run())
        #self.measurements = self.subsetting(self.measurements)
        self.output()

    def subsetting(self, cursor):
        '''
        filter through documents in vdb to return subsets of sequence
        '''
        result = cursor
        print("Documents in table before subsetting: " + str(len(result)))
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

    def write_txt(self, measurements, fname):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for meas in self.measurements:
                for field in self.txt_fields:
                    if field in meas and meas[field] is not None:
                        handle.write(str(meas[field]) + '\t')
                    else:
                        handle.write('?' + '\t')
                handle.write('\n')
            handle.close()
            print("Wrote to " + fname)

    def output(self):
        if self.ftype == 'json':
            self.write_json(self.measurements, self.path+self.fname)
        else:
            self.write_txt(self.measurements, self.path+self.fname)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = vdb_download(**args.__dict__)
    connVDB.download()