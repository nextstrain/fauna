import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
import argparse
from tdb_parse import tdb_parse

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='tdb', help="database to upload to")
parser.add_argument('-v', '--virus', default='H3N2', help="virus table to interact with, ie \'H3N2\' or \'Yam\' for Flu")
parser.add_argument('--path', default=None, help="path to fasta file, default is \"data/virus/\"")
parser.add_argument('--overwrite', default=False, action="store_true",  help ="Overwrite fields that are not none")
parser.add_argument('--exclusive', default=True, action="store_false",  help ="download all docs in db to check before upload")
parser.add_argument('--upload', default=False, action="store_true",  help ="If included, actually upload documents, otherwise test parsing")
parser.add_argument('--host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class tdb_upload(tdb_parse):
    def __init__(self, **kwargs):
        tdb_parse.__init__(self, **kwargs)
        if 'database' in kwargs:
            self.database = kwargs['database']
        if 'virus' in kwargs:
            self.virus = kwargs['virus'].lower()
        if 'overwrite' in kwargs:
            self.overwrite = kwargs['overwrite']
        if 'exclusive' in kwargs:
            self.exclusive = kwargs['exclusive']
        if 'upload' in kwargs:
            self.upload_docs = kwargs['upload']
        if 'path' in kwargs:
            self.path = kwargs['path']
        if self.path is None:
            self.path = "tdb/data/" + self.virus + "/"
        if not os.path.isdir(self.path):
            os.makedirs(self.path)

        # fields that are needed to upload
        self.upload_fields = ['virus', 'serum', 'titer', 'date_modified', 'source', 'ferret_id', 'passage'] #index too but assign after checking
        self.optional_fields = ['date', 'ref']
        self.overwritable_fields = ['titer', 'date', 'ref']
        self.index_fields = ['virus', 'serum', 'ferret_id', 'source', 'passage']

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
        self.test_virus_strains = set()
        self.indexes = set()
        self.passage = set()
        self.different_date_format = ['NIMR-REPORT-FEB2010_03.CSV', 'NIMR-REPORT-FEB2010_06.CSV', 'NIMR-REPORT-FEB2010_05.CSV',
                                      'NIMR-REPORT-SEP2009_03.CSV', 'NIMR-REPORT-FEB2010_04.CSV', 'NIMR_FEB2010_15.CSV', 'NIMR_FEB2010_16.CSV']
        self.new_different_date_format = set()

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
            raise Exception("No table exists yet for " + str(self.virus))

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
        self.create_index()
        print('Total number of indexes', len(self.indexes), 'Total number of measurements', len(self.measurements))
        if self.upload_docs:
            self.upload_documents()


    def format(self):
        '''
        format virus information in preparation to upload to database table
        '''
        print('Formatting for upload')
        for meas in self.measurements:
            self.format_date(meas)
            #self.assign_numdate(meas)
            self.format_passage(meas)
            self.format_id(meas)
            self.format_ref(meas)
            self.format_titer(meas)
            if meas['ref'] == True:
                self.ref_serum_strains.add(meas['serum'])
                self.ref_virus_strains.add(meas['virus'])
            if meas['ref'] == False:
                self.test_virus_strains.add(meas['virus'])
        if len(self.new_different_date_format) > 0:
            print("Found files that had a different date format, need to add to self.different_date_format")
            print(self.new_different_date_format)
        self.check_strain_names()

    def check_strain_names(self):
        '''
        Check that strain names for serum have been parsed correctly, filter out those that haven't
        Virus strain names are assumed to be parsed correctly
        '''
        weird_serum = []
        if self.ref_serum_strains != self.ref_virus_strains:
            for meas in self.measurements:
                item = meas['serum'], meas['source']
                if meas['serum'] in (self.ref_serum_strains - self.ref_virus_strains - self.test_virus_strains) and item not in weird_serum:
                    weird_serum.append(item)
        if len(weird_serum) > 0:
            print("Warning: different serum and virus strain names, check for missing abbreviations to add to tdb_parse")
            print("Will filter out measurements with different serum names")
            print("Actual reference virus names")
            print(self.ref_virus_strains)
            print("Odd reference serum names")
            print(self.ref_serum_strains - self.ref_virus_strains - self.test_virus_strains)
            print("Odd reference serum name locations")
            print(weird_serum)
            #print("Test reference virus names")
            #print(self.test_virus_strains)

    def format_date(self, meas):
        '''
        Format date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        Input date could be YYYY_MM_DD, reformat to YYYY-MM-DD
        '''
        # ex. 2002_04_25 to 2002-04-25
        lookup_month = {'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04', 'May': '05', 'Jun': '06',
                        'Jul': '07', 'Aug': '08', 'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
        try:
            meas['date'] = re.sub(r'_', r'-', meas['date']).strip()
        except:
            meas['date'] = ''
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
        elif re.match(r'^\d+/\d+/\d\d$', meas['date']):
            try:
                if meas['source'] not in self.different_date_format:
                    date = datetime.datetime.strptime(meas['date'], '%m/%d/%y').date()
                else:
                    date = datetime.datetime.strptime(meas['date'], '%d/%m/%y').date()
                meas['date'] = date.strftime('%Y-%m-%d')
            except:
                meas['date'] = None
                self.new_different_date_format.add(meas['source'])
        elif re.match(r'^\d+/\d+/\d\d\d\d$', meas['date']):
            try:
                if meas['source'] not in self.different_date_format:
                    date = datetime.datetime.strptime(meas['date'], '%m/%d/%Y').date()
                else:
                    date = datetime.datetime.strptime(meas['date'], '%d/%m/%Y').date()
                meas['date'] = date.strftime('%Y-%m-%d')
            except:
                meas['date'] = None
                self.new_different_date_format.add(meas['source'])
        elif re.match(r'[a-zA-Z][a-zA-Z][a-zA-Z]\s\d\d\d\d', meas['date']):  #Jan 2012
            try:
                date = datetime.datetime.strptime(meas['date'], '%b %Y').date()
                meas['date'] = date.strftime('%Y-%m') + '-XX'
            except:
                print("Couldn't parse as datetime object", meas['date'], meas['source'])
                meas['date'] = None
        elif re.match(r'(\d+)-[a-zA-Z][a-zA-Z][a-zA-Z]', meas['date']):  #12-Jan, 9-Jun
            try:
                year = re.match(r'(\d+)-[a-zA-Z][a-zA-Z][a-zA-Z]', meas['date']).group(1)
                if len(str(year)) == 1:
                    meas['date'] = '0' + meas['date']
                date = datetime.datetime.strptime(meas['date'], '%y-%b').date()
                meas['date'] = date.strftime('%Y-%m') + '-XX'
                print(meas['date'])
            except:
                print("Couldn't parse as datetime object", meas['date'], meas['source'])
                meas['date'] = None
        elif re.match(r'[a-zA-Z][a-zA-Z][a-zA-Z]-(\d+)', meas['date']):  #Nov-09, Jan-2012
            try:
                year = re.match(r'[a-zA-Z][a-zA-Z][a-zA-Z]-(\d+)', meas['date']).group(1)
                if len(str(year)) == 4:
                    date = datetime.datetime.strptime(meas['date'], '%b-%Y').date()
                else:
                    date = datetime.datetime.strptime(meas['date'], '%b-%y').date()
                meas['date'] = date.strftime('%Y-%m') + '-XX'
                print(meas['date'])
            except:
                print("Couldn't parse as datetime object", meas['date'], meas['source'])
                meas['date'] = None
        elif meas['date'].lower() == 'unknown' or meas['date'].lower() == 'd/m unknown':
            meas['date'] = None
        elif meas['date'] == '':
            meas['date'] = None
        elif meas['date'] == 'Si':
            meas['date'] = '2009-07-01'
            #print(meas['virus'], 'A/BAYERN/69/2009', 'date assigned to 2009-07-01')
        else:
            print("Couldn't reformat this date: \'" + meas['date'] + "\'", meas['source'], meas['serum'], meas['virus'])
            meas['date'] = None

    def assign_numdate(self, meas, format = '%Y-%m-%d'):
        '''
        Takes a calendar date and a numerical dates in terms of years
        If date is missing precision, still return numerical date
        '''
        meas['numdate'] = meas['date']
        if meas['date'] is not None:
            date = meas['date'].replace("XX", "01")
            if isinstance(date, basestring):
                try:
                    date = datetime.datetime.strptime(date, format).date()
                    meas['numdate'] = round(date.year + 1.0*date.timetuple().tm_yday/365.25, 3)
                except:
                    print("Couldn't parse date, may need to fix manually", meas['date'], meas['source'], meas['serum'], meas['virus'])

    def format_passage(self, meas):
        '''
        Format passage attribute
        '''

    def format_id(self, meas):
        '''
        Format ferret id attribute
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

    def format_titer(self, meas):
        '''
        Format titer number attribute
        '''

    def create_index(self):
        '''
        create unique key for storage in rethinkdb
        '''
        for meas in self.measurements:
            index = [meas[field] for field in self.index_fields]
            meas['index'] = index
            if str(index) in self.indexes:
                print("Repeat Index: " + str(index) + str(meas['titer']))
            else:
                self.indexes.add(str(index))

    def filter(self):
        '''
        Filter out viruses without correct dating format or without region specified
        Check  optional and upload attributes
        '''
        print(len(self.measurements), " measurements before filtering")
        self.check_optional_attributes()
        self.measurements = filter(lambda meas: isinstance(meas['ref'], (bool)), self.measurements)
        self.measurements = filter(lambda meas: meas['serum'] in self.ref_virus_strains or meas['serum'] in self.test_virus_strains, self.measurements)
        self.measurements = filter(lambda meas: self.check_upload_attributes(meas), self.measurements)
        print(len(self.measurements), " measurements after filtering")

    def check_optional_attributes(self):
        '''
        Create and assign 'None' to optional attributes that don't exist
        '''
        for meas in self.measurements:
            for atr in self.optional_fields + self.upload_fields:
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
            #print("Missing required attributes:" + str(missing_attributes) + ", measurement: " + meas['virus'], meas['serum'], meas['source'])
            return False
        else:
            return True

    def upload_documents(self):
        '''
        Insert viruses into collection
        '''
        if self.exclusive:
            db_measurements = list(r.db(self.database).table(self.virus).run())
            uploaded_indexes = {" ".join([str(i) for i in db_meas['index']]): db_meas for db_meas in db_measurements}
            check_update_measurements = []
            upload_measurements = []
            added_to_upload_measurements = set()
            for meas in self.measurements:
                index_test = " ".join([str(i) for i in meas['index']])
                if index_test not in uploaded_indexes.keys() and index_test not in added_to_upload_measurements:
                    upload_measurements.append(meas)
                    added_to_upload_measurements.add(index_test)
                elif index_test in uploaded_indexes.keys():
                    check_update_measurements.append(meas)
            print("Inserting ", len(upload_measurements), "measurements into database", self.virus)
            r.table(self.virus).insert(upload_measurements).run()
            print("Checking for updates to ", len(check_update_measurements), "measurements in database", self.virus)
            for meas in check_update_measurements:
                pass
                self.updated = False
                index_test = " ".join([str(i) for i in meas['index']])
                self.update_document_meta(uploaded_indexes[index_test], meas)
        else:
            print("Uploading " + str(len(self.measurements)) + " measurements to the table")
            for meas in self.measurements:
                try:
                    document = r.table(self.virus).get(meas['index']).run()
                except:
                    print(meas)
                    raise Exception("Couldn't retrieve this measurement")
                # Virus doesn't exist in table yet so add it
                if document is None:
                    #print("Inserting " + meas['index'] + " into database")
                    try:
                        r.table(self.virus).insert(meas).run()
                    except:
                        print(meas)
                        raise Exception("Couldn't insert this measurement")
                # Virus exists in table so just add sequence information and update meta data if needed
                else:
                    self.updated = False
                    self.update_document_meta(document, meas)

    def update_document_meta(self, document, meas):
        '''
        if overwrite is false update doc_virus information only if virus info is different and not null
        if overwrite is true update doc_virus information if virus info is different
        '''
        for field in self.overwritable_fields:
            # update if overwrite and anything
            # update if !overwrite only if document[field] is not none
            if field not in document:
                if field in meas:
                    print("Creating measurement field ", field, " assigned to ", meas[field])
                    r.table(self.virus).get(meas['index']).update({field: meas[field]}).run()
                    document[field] = meas[field]
                    self.updated = True
            elif (self.overwrite and document[field] != meas[field]) or (not self.overwrite and document[field] is None and document[field] != meas[field]):
                if field in meas:
                    print("Updating measurement field " + str(field) + ", from \"" + str(document[field]) + "\" to \"" + meas[field]) + "\""
                    r.table(self.virus).get(meas['index']).update({field: meas[field]}).run()
                    document[field] = meas[field]
                    self.updated = True
        if self.updated:
            r.table(self.virus).get(meas['index']).update({'date_modified': meas['date_modified']}).run()
            document['date_modified'] = meas['date_modified']

if __name__=="__main__":
    args = parser.parse_args()
    connTDB = tdb_upload(**args.__dict__)
    connTDB.upload()