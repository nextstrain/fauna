import os, re, time, datetime, csv, sys, json
import rethinkdb as r
from Bio import SeqIO
import argparse
from parse import parse
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='tdb', help="database to upload to")
parser.add_argument('-v', '--virus', default='flu', help="virus table to interact with, ie Flu")
parser.add_argument('--subtype', default=None, help="subtype of virus, ie h1n1pdm, vic, yam, h3n2")
parser.add_argument('--host', default='human', help="host of virus, ie human, swine")
parser.add_argument('--path', default=None, help="path to fasta file, default is \"data/virus/\"")
parser.add_argument('--overwrite', default=False, action="store_true",  help ="Overwrite fields that are not none")
parser.add_argument('--exclusive', default=True, action="store_false",  help ="download all docs in db to check before upload")
parser.add_argument('--replace', default=False, action="store_true",  help ="If included, delete all documents in table")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('--preview', default=False, action="store_true",  help ="If included, preview a virus document to be uploaded")

class upload(parse, flu_upload):
    def __init__(self, database, virus, subtype, **kwargs):
        parse.__init__(self, **kwargs)
        flu_upload.__init__(self, database=database, virus=virus)
        self.virus = virus.lower()
        self.subtype = subtype.lower()
        self.database = database.lower()
        self.uploadable_databases = ['tdb', 'test_tdb', 'test']
        if self.database not in self.uploadable_databases:
            raise Exception("Cant upload to this database: " + self.database, "add to list of databases allowed", self.uploadable_databases)
        self.flu_upload = flu_upload(database=self.database, virus=self.virus)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.virus)
        self.strain_fix_fname = "source-data/flu_strain_name_fix.tsv"
        self.HI_strain_fix_fname = "source-data/HI_flu_strain_name_fix.tsv"
        self.HI_ref_name_abbrev_fname = "source-data/HI_ref_name_abbreviations.tsv"

        # fields that are needed to upload
        self.upload_fields = ['virus_strain', 'serum_strain', 'titer', 'timestamp', 'source', 'ferret_id', 'passage', 'subtype',
                              'host'] #index too but assign after checking
        self.optional_fields = ['date', 'ref']
        self.overwritable_fields = ['titer', 'date', 'ref']
        self.index_fields = ['virus_strain', 'serum_strain', 'ferret_id', 'source', 'passage', 'subtype', 'host']
        self.ref_virus_strains = set()
        self.ref_serum_strains = set()
        self.test_virus_strains = set()
        self.indexes = set()
        self.passage = set()
        self.strain_names = set()
        self.different_date_format = ['NIMR-REPORT-FEB2010_03.CSV', 'NIMR-REPORT-FEB2010_06.CSV', 'NIMR-REPORT-FEB2010_05.CSV', 'NIMR_Feb2010_15.csv',
                                      'NIMR-REPORT-SEP2009_03.CSV', 'NIMR-REPORT-FEB2010_04.CSV', 'NIMR_FEB2010_15.CSV', 'NIMR_FEB2010_16.CSV', 'NIMR_Feb2010_16.csv',
                                      'NIMR-report-Feb2010_03.csv', 'NIMR-report-Feb2010_06.csv', 'NIMR-report-Feb2010_05.csv', 'NIMR-report-Sep2009_03.csv', 'NIMR-report-Feb2010_04.csv']
        self.new_different_date_format = set()
        self.fix = set()

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to TDB")
        measurements = self.parse(**kwargs)
        print('Formatting documents for upload')
        self.format_measurements(measurements, **kwargs)
        measurements = self.filter(measurements)
        measurements = self.create_index(measurements)
        print('Total number of indexes', len(self.indexes), 'Total number of measurements', len(measurements))
        if not preview:
            self.upload_documents(measurements, **kwargs)
        else:
            print("Titer Measurements:")
            print(json.dumps(measurements[0], indent=1))
            print("Remove \"--preview\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def format_measurements(self, measurements, vdb_database='test_vdb', vdb_table='flu_viruses', **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        self.fix_whole_name.update(self.define_strain_fixes(self.HI_strain_fix_fname))
        self.HI_ref_name_abbrev =self.define_strain_fixes(self.HI_ref_name_abbrev_fname)
        self.define_location_fixes("source-data/flu_fix_location_label.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        print("Using " + vdb_database + "." + vdb_table + " to adjust strain names to sequence strains")
        self.vdb_strains = self.relaxed_keys(set(list(r.db(vdb_database).table(vdb_table).get_field('strain').run())))
        self.unknown_strains = 0
        self.adjusted_strains = 0
        self.adjusted = set()
        self.known_strains = 0
        for meas in measurements:
            meas['virus_strain'], meas['original_virus_strain'] = self.fix_name(self.HI_fix_name(meas['virus_strain'], serum=False))
            self.correct_strain_format(meas['virus_strain'], meas['original_virus_strain'])
            meas['serum_strain'], meas['original_serum_strain'] = self.fix_name(self.HI_fix_name(meas['serum_strain'], serum=True))
            self.correct_strain_format(meas['serum_strain'], meas['original_serum_strain'])
            relax_strain = self.relax_name(meas['virus_strain'])
            if relax_strain not in self.vdb_strains:  # no sequence for this strain
                self.unknown_strains += 1
            elif self.vdb_strains[relax_strain] != meas['virus_strain']:  # sequence, but strain name needs to be adjusted
                self.adjusted.add(meas['virus_strain'] + " : " + self.vdb_strains[relax_strain])
                meas['virus_strain'] = self.vdb_strains[relax_strain]
                self.adjusted_strains += 1
            else:  # sequence known for this strain name
                self.known_strains += 1
            self.test_location(meas['virus_strain'])
            self.test_location(meas['serum_strain'])
            self.add_attributes(meas, **kwargs)
            self.format_date(meas)
            self.format_passage(meas)
            self.format_id(meas)
            self.format_ref(meas)
            self.format_titer(meas)
            if meas['ref'] == True:
                self.ref_serum_strains.add(meas['serum_strain'])
                self.ref_virus_strains.add(meas['virus_strain'])
            if meas['ref'] == False:
                self.test_virus_strains.add(meas['virus_strain'])
            self.rethink_io.check_optional_attributes(meas, self.optional_fields)
        if len(self.new_different_date_format) > 0:
            print("Found files that had a different date format, need to add to self.different_date_format")
            print(self.new_different_date_format)
        self.check_strain_names(measurements)
        print("No sequence for these virus strains", self.unknown_strains)
        print("Strain name matched vdb strain initially", self.known_strains)
        print("Adjusted strain names to match vdb strain", self.adjusted_strains)
        print("Titer measurement strain name : Sequence strain name")
        for item in sorted(self.adjusted):
            print(item)
        return measurements

    def HI_fix_name(self, name, serum):
        lookup_month = {'Jan': '1', 'Feb': '2', 'Mar': '3', 'Apr': '4', 'May': '5', 'Jun': '6',
                            'Jul': '7', 'Aug': '8', 'Sep': '9', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
        name = name.replace('H1N1', '').replace('H5N6', '').replace('H3N2', '').replace('Human', '')\
            .replace('human', '').replace('//', '/').replace('.', '').replace(',', '').replace('&', '').replace(' ', '')\
            .replace('\'', '').replace('>', '').replace('-like', '').replace('+', '').replace('*', '')
        if re.match('[0-9]{1,2}([A|B]/.+)', name):  # 12B/Estonia/55669/2011 -> B/Estonia/55669/2011, 2B/Estonia/55669/2011 -> B/Estonia/55669/2011
            name = re.match('[0-9]{1,2}([A|B]/.+)', name).group(1)
        if serum and re.match(r'([A|B]/[A-Za-z-]+)[0-9]+(/[A-Za-z0-9_-]+/[0-9]{2,4})', name):  # B/Bris13/60/08 -> B/Bris/60/08, B/Fl13/6-Apr -> B/Fl/6-Apr
            name = re.match(r'([A|B]/[A-Za-z-]+)[0-9]+(/[A-Za-z0-9_-]+/[0-9]{2,4})', name).group(1) + re.match(r'([A|B]/[A-Za-z-]+)[0-9]+(/[A-Za-z0-9_-]+/[0-9]{2,4})', name).group(2)
        if serum and re.match(r'([A|B]/[A-Za-z-]+)[0-9]+(/[0-9]{1,2}-[A-Za-z]{3})', name):  # B/Fl1/6-Apr -> B/Fl/6-Apr, B/Bris1/7-Mar -> B/Bris/7-Mar
            name = re.match(r'([A|B]/[A-Za-z-]+)[0-9]+(/[0-9]{1,2}-[A-Za-z]{3})', name).group(1) + re.match(r'([A|B]/[A-Za-z-]+)[0-9]+(/[0-9]{1,2}-[A-Za-z]{3})', name).group(2)
        if re.match(r'([A|B]/[A-Za-z-]+/)([0-9]{1,2}-[A-Za-z]{3})', name):  # B/Fl/6-Apr -> B/Fl/4/2006, B/Stock/11-Dec -> B/Stock/12/2011
            try:
                date_pattern = re.match(r'([A|B]/[A-Za-z-]+/)([0-9]{1,2}-[A-Za-z]{3})', name).group(2).split('-')
                year = int(date_pattern[0])
                if year>68:
                    year=1900+year
                else:
                    year=2000+year
                date = lookup_month[date_pattern[1]] + "/" + str(year)
                name = re.match(r'([A|B]/[A-Za-z-]+/)([0-9]{1,2}-[A-Za-z]{3})', name).group(1) + date
            except:
                pass
        if re.match(r'([A|B]/[A-Za-z-]+/)([A-Za-z]{3}-[0-9]{1,2})', name):  # B/SHANDONG/JUL-97 -> B/SHANDONG/7/1997, A/NewJersey/8/1976
            try:
                date_pattern = re.match(r'([A|B]/[A-Za-z-]+/)([A-Za-z]{3}-[0-9]{1,2})', name).group(2).split('-')
                year = int(date_pattern[1])
                if year>68:
                    year=1900+year
                else:
                    year=2000+year
                date = lookup_month[date_pattern[0]] + "/" + str(year)
                name = re.match(r'([A|B]/[A-Za-z-]+/)([A-Za-z]{3}-[0-9]{1,2})', name).group(1) + date
            except:
                pass
        if re.match(r'([A|B]/)([A-Za-z-]+)(/.+)', name):
            match = re.match(r'([A|B]/)([A-Za-z-]+)(/.+)', name)
            abbrev = match.group(2).upper()
            if abbrev in self.HI_ref_name_abbrev:
                name = match.group(1) + self.HI_ref_name_abbrev[abbrev] + match.group(3)
        return name

    def test_location(self, strain):
        if isinstance(strain, basestring) and "/" in strain:
            location = strain.split('/')[1]
            if self.determine_location(location) is None:
                print("Couldn't determine location for this strain, consider adding to flu_fix_location_label.tsv", location, strain)

    def add_attributes(self, meas, host, **kwargs):
        '''
        Add attributes to virus
        '''
        meas['virus'] = self.virus.lower()
        meas['subtype'] = self.subtype.lower()
        meas['host'] = host.lower()
        meas['timestamp'] = self.rethink_io.get_upload_timestamp()

    def check_strain_names(self, measurements):
        '''
        Check that strain names for serum have been parsed correctly, filter out those that haven't
        Virus strain names are assumed to be parsed correctly
        '''
        weird_serum = []
        if self.ref_serum_strains != self.ref_virus_strains:
            for meas in measurements:
                item = meas['serum_strain'], meas['source']
                if meas['serum_strain'] in (self.ref_serum_strains - self.ref_virus_strains - self.test_virus_strains) and item not in weird_serum:
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
            except:
                print("Couldn't parse as datetime object", meas['date'], meas['source'])
                meas['date'] = None
        elif meas['date'].lower() == 'unknown' or meas['date'].lower() == 'd/m unknown':
            meas['date'] = None
        elif meas['date'] == '':
            meas['date'] = None
        else:
            print("Couldn't reformat this date: \'" + meas['date'] + "\'", meas['source'], meas['serum_strain'], meas['virus_strain'])
            meas['date'] = None

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
        if 'ref' in meas:
            if meas['ref'].lower() == 'ref':
                meas['ref'] = True
            elif meas['ref'].lower() == 'test':
                meas['ref'] = False
            else:
                meas['ref'] = None
                print("Couldn't parse reference status", meas['ref'])
        else:
            meas['ref'] = None

    def format_titer(self, meas):
        '''
        Format titer number attribute
        '''

    def create_index(self,  measurements, output=False):
        '''
        create unique key for storage in rethinkdb
        '''
        for meas in measurements:
            index = [meas[field] for field in self.index_fields]
            meas['index'] = index
            if str(index) in self.indexes:
                if output:
                    print("Repeat Index: " + str(index) + str(meas['titer']))
            else:
                self.indexes.add(str(index))
        return measurements

    def filter(self, measurements):
        '''
        Filter out viruses without correct dating format or without region specified
        Check  optional and upload attributes
        '''
        print(len(measurements), " measurements before filtering")
        print("Filtering out measurements whose serum strain is not paired with a ref or test virus strain ensuring proper formatting")
        measurements = filter(lambda meas: meas['serum_strain'] in self.ref_virus_strains or meas['serum_strain'] in self.test_virus_strains, measurements)
        print("Filtering out measurements missing required fields")
        measurements = filter(lambda meas: self.rethink_io.check_required_attributes(meas, self.upload_fields, self.index_fields), measurements)
        print(len(measurements), " measurements after filtering")
        return measurements

    def upload_documents(self, measurements, replace=False, exclusive=True, **kwargs):
        '''
        Insert viruses into collection
        '''
        if replace:
            print("Deleting documents in database:", self.database, "table:", self.virus)
            r.table(self.virus).delete().run()
        if exclusive:
            db_measurements = list(r.db(self.database).table(self.virus).run())
            uploaded_indexes = {" ".join([str(i) for i in db_meas['index']]): db_meas for db_meas in db_measurements}
            check_update_measurements = []
            upload_measurements = []
            added_to_upload_measurements = set()
            for meas in measurements:
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
                self.updated = False
                index_test = " ".join([str(i) for i in meas['index']])
                self.update_document_meta(uploaded_indexes[index_test], meas)
        else:
            print("Uploading " + str(len(measurements)) + " measurements to the table")
            for meas in measurements:
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
                    self.update_document_meta(document, meas, **kwargs)

    def upload_to_rethinkdb(self, database, table, documents, conflict_resolution, **kwargs):
        optimal_upload = 200
        if len(documents) > optimal_upload:
            list_documents = [documents[x:x+optimal_upload] for x in range(0, len(documents), optimal_upload)]
        else:
            list_documents = [documents]
        print("Uploading to rethinkdb in " + str(len(list_documents)) + " batches")
        for list_docs in list_documents:
            try:
                r.table(table).insert(list_docs, conflict=conflict_resolution).run()
            except:
                raise Exception("Couldn't insert new documents into database", database + "." + table)

    def update_document_meta(self, document, meas, overwrite=False, **kwargs):
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
            elif (overwrite and document[field] != meas[field]) or (not overwrite and document[field] is None and document[field] != meas[field]):
                if field in meas:
                    print("Updating measurement field " + str(field) + ", from \"" + str(document[field]) + "\" to \"" + meas[field]) + "\""
                    r.table(self.virus).get(meas['index']).update({field: meas[field]}).run()
                    document[field] = meas[field]
                    self.updated = True
        if self.updated:
            r.table(self.virus).get(meas['index']).update({'timestamp': meas['timestamp']}).run()
            document['timestamp'] = meas['timestamp']

    def relaxed_keys(self, indexes):
        '''
        Create dictionary from relaxed vdb strain names to actual vdb strain names.
        '''
        strains = {}
        for index in indexes:
            strains[self.relax_name(index)] = index
        return strains

    def relax_name(self, name):
        '''
        Return the relaxed strain name to compare with
        '''
        split_name = name.split('/')
        for index, split in enumerate(split_name):  # A/Guangdong-Yuexiu/SWL1651/2013 -> A/Guangdong-Yuexiu/1651/2013
            if index >= 2:
                split_name[index] = filter(lambda x: x.isdigit(), split)
        name = "/".join(split_name)
        name = re.sub(r"-", '', name)
        name = re.sub(r"_", '', name)
        name = re.sub(r"/", '', name)
        return name

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "tdb/data/" + args.subtype + "/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = upload(**args.__dict__)
    connTDB.upload(**args.__dict__)