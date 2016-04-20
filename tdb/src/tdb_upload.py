import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
import argparse
from tdb_parse import tdb_parse

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='tdb', help="database to upload to")
parser.add_argument('-v', '--virus', default='H3N2', help="virus table to interact with, ie \'H3N2\' or \'Yam\' for Flu")
parser.add_argument('--fname', help="input file name")
parser.add_argument('--path', default=None, help="path to fasta file, default is \"data/virus/\"")
parser.add_argument('--host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class tdb_upload(tdb_parse):
    def __init__(self, **kwargs):
        tdb_parse.__init__(self, **kwargs)

        if 'database' in kwargs:
            self.database = kwargs['database']
        if 'virus' in kwargs:
            self.virus = kwargs['virus'].upper()
        if 'fname' in kwargs:
            self.fname = kwargs['fname']

        if 'path' in kwargs:
            self.path = kwargs['path']
        if self.path is None:
            self.path = "tdb/data/" + self.virus + "/"
        if not os.path.isdir(self.path):
            os.makedirs(self.path)

        # fields that are needed to upload
        self.upload_fields = ['virus', 'serum', 'titer', 'ferret_id', 'date_modified', 'source']
        self.optional_fields = ['date', 'passage', 'group', 'ref']

        if 'host' in kwargs:
            self.host = kwargs['host']
        if 'RETHINK_HOST' in os.environ and self.host is None:
            self.host = os.environ['RETHINK_HOST']
        if self.host is None:
            raise Exception("Missing rethink host")

        if 'auth_key' in kwargs:
            self.auth_key = kwargs['auth_key']
        if 'RETHINK_AUTH_KEY' in os.environ and self.auth_key is None:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")

        self.connect_rethink()
        self.ref_virus_strains = set()
        self.ref_serum_strains = set()

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

    def get_upload_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))

    def upload(self):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to TDB")
        self.parse()
        self.format()
        self.filter()
        self.upload_documents()

    def format(self):
        '''
        format virus information in preparation to upload to database table
        '''
        print('Formatting for upload')
        for meas in self.measurements:
            self.format_date(meas)
            self.format_passage(meas)
            self.format_id(meas)
            self.format_group(meas)
            self.format_ref(meas)
            self.create_index(meas)
            if meas['ref'] == True:
                self.ref_serum_strains.add(meas['serum'])
                self.ref_virus_strains.add(meas['virus'])
        self.check_strain_names()

    def check_strain_names(self):
        '''
        Check that strain names for serum have been parsed correctly, filter out those that haven't
        Virus strain names are assumed to be parsed correctly
        '''
        if self.ref_serum_strains != self.ref_virus_strains:
            print("Warning: different serum and virus strain names, check for missing abbreviations to add to tdb_parse")
            print("Will filter out measurements with different serum names")
        print("Odd reference serum names")
        print(self.ref_serum_strains - self.ref_virus_strains)
        print("Actual reference virus names")
        print(self.ref_virus_strains)

    def format_date(self, meas):
        '''
        Format date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        Input date could be YYYY_MM_DD, reformat to YYYY-MM-DD
        '''
        # ex. 2002_04_25 to 2002-04-25
        lookup_month = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06',
                        'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
        meas['date'] = re.sub(r'_', r'-', meas['date']).strip()
        # ex. 2002 (Month and day unknown)
        if re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', meas['date']):
            pass
        elif re.match(r'\d\d\d\d\s\(Month\sand\sday\sunknown\)', meas['date']):
            meas['date'] = meas['date'][0:4] + "-XX-XX"
        # ex. 2009-06 (Day unknown)
        elif re.match(r'\d\d\d\d-\d\d\s\(Day\sunknown\)', meas['date']):
            meas['date'] = meas['date'][0:7] + "-XX"
        elif re.match(r'\d\d\d\d-\d\d', meas['date']):
            meas['date'] = meas['date'][0:7] + "-XX"
        elif re.match(r'\d\d\d\d', meas['date']):
            meas['date'] = meas['date'][0:4] + "-XX-XX"
        elif re.match(r'\d\d/\d\d/\d\d\d\d', meas['date']):
            meas['date'] = meas['date'][6:10] + "-" + meas['date'][3:5] + "-" + meas['date'][0:2]
        elif re.match(r'[a-zA-Z][a-zA-Z][a-zA-Z]-\d\d', meas['date']):
            meas['date'] = '20' + meas['date'][4:6] + '-' + lookup_month[meas['date'][0:3]] + "-XX"
        elif meas['date'] == 'unknown':
            meas['date'] = None
        elif meas['date'] == '':
            meas['date'] = None
        elif meas['date'] == 'Si':
            meas['date'] = '2009-07-01'
            print(meas['virus'], 'A/BAYERN/69/2009', 'date assigned to 2009-07-01')
        else:
            print("Couldn't reformat this date: \'" + meas['date'] + "\'")

    def format_passage(self, meas):
        '''
        Format passage attribute
        '''

    def format_id(self, meas):
        '''
        Format ferret id attribute
        '''

    def format_group(self, meas):
        '''
        Format genetic group attribute
        '''

    def format_ref(self, meas):
        '''
        Format ref/test attribute
        '''
        if meas['ref'] == 'REF':
            meas['ref'] = True
        elif meas['ref'] == 'TEST':
            meas['ref'] = False
        else:
            print("Couldn't parse reference status of this measurement", meas)

    def create_index(self, meas):
        '''
        create unique key for storage in rethinkdb
        '''
        index = [meas['virus'], meas['serum'], meas['ferret_id']]
        meas['index'] = index

    def filter(self):
        '''
        Filter out viruses without correct dating format or without region specified
        Check  optional and upload attributes
        '''
        print(len(self.measurements), " measurements before filtering")
        self.check_optional_attributes()
        self.measurements = filter(lambda meas: isinstance(meas['ref'], (bool)), self.measurements)
        self.measurements = filter(lambda meas: meas['serum'] in self.ref_virus_strains, self.measurements)
        self.measurements = filter(lambda meas: self.check_upload_attributes(meas), self.measurements)
        print(len(self.measurements), " measurements after filtering")

    def check_optional_attributes(self):
        '''
        Create and assign 'None' to optional attributes that don't exist
        '''
        for meas in self.measurements:
            for atr in self.optional_fields:
                if atr not in meas:
                    meas[atr] = None
                elif meas[atr] == '':
                    meas[atr] = None

    def check_upload_attributes(self, meas):
        '''
        Checks that required upload attributes are present and not equal to None for given virus
        :return: returns true if it has all required upload attributes, else returns false and prints missing attributes
        '''
        missing_attributes = []
        for atr in self.upload_fields:
            if atr not in meas or meas[atr] is None:
                missing_attributes.append(atr)
        if len(missing_attributes) > 0:
            print("This measurement is missing a required attribute and will be removed from upload sequences")
            print(meas['virus'], meas['serum']), meas['id']
            print("Missing attributes: " + str(missing_attributes))
            return False
        else:
            return True

    def upload_documents(self):
        '''
        Insert viruses into collection
        '''
        print("Uploading " + str(len(self.measurements)) + " measurements to the table")
        for meas in self.measurements:
            document = r.table(self.virus).get(meas['index']).run()
            # Virus doesn't exist in table yet so add it
            if document is None:
                print("Inserting " + meas['index'] + " into database")
                r.table(self.virus).insert(meas).run()
            # Virus exists in table so just add sequence information and update meta data if needed
            else:
                pass
                #self.updated = False
                #self.update_document_sequence(document, meas)
                #self.update_document_meta(document, meas)

if __name__=="__main__":
    args = parser.parse_args()
    run = tdb_upload(**args.__dict__)
    run.upload()