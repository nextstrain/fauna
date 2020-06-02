from __future__ import print_function
import os, re, time, datetime, csv, sys, json
import rethinkdb as r
import hashlib
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
parser.add_argument('--path', default='data/', help="path to input file, default is \"data/\"")
parser.add_argument('--fstem', help="input file stem")
parser.add_argument('--ftype', default='flat', help="input file format, default \"flat\", other is \"tables\"")
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
        self.table = self.virus
        if subtype is not None:
            self.subtype = subtype.lower()
        self.database = database.lower()
        self.uploadable_databases = ['tdb', 'test_tdb', 'test_tdb_2', 'cdc_tdb', 'test_cdc', 'vidrl_tdb', 'test_vidrl', 'crick_tdb', 'test_crick', 'niid_tdb', 'test_niid']
        if self.database not in self.uploadable_databases:
            raise Exception("Can't upload to this database: " + self.database, "add to list of databases allowed", self.uploadable_databases)
        self.flu_upload = flu_upload(database=self.database, virus=self.virus)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.virus)
        self.strain_fix_fname = "source-data/flu_strain_name_fix.tsv"
        self.HI_ref_name_abbrev_fname = "source-data/HI_ref_name_abbreviations.tsv"

        # fields that are needed to upload
        self.upload_fields = ['virus_strain', 'serum_strain', 'serum_id', 'assay_type', 'assay_date', 'titer', 'timestamp']
        self.optional_fields = ['date', 'ref', 'virus_cdc_id', 'virus_passage', 'virus_passage_category', 'subtype', 'serum_host', 'serum_passage', 'serum_passage_category']
        self.overwritable_fields = ['titer', 'date', 'ref']
        self.index_fields = ['virus_strain', 'serum_strain', 'serum_id', 'source', 'virus_passage_category', 'serum_passage_category', 'assay_type', 'assay_date']
        self.ref_virus_strains = set()
        self.ref_serum_strains = set()
        self.test_virus_strains = set()
        self.indexes = set()
        self.virus_passage = set()
        self.serum_passage = set()
        self.strain_names = set()
        self.assay_type = set()
        self.different_date_format = ['NIMR-REPORT-FEB2010_03.CSV', 'NIMR-REPORT-FEB2010_06.CSV', 'NIMR-REPORT-FEB2010_05.CSV', 'NIMR_Feb2010_15.csv',
                                      'NIMR-REPORT-SEP2009_03.CSV', 'NIMR-REPORT-FEB2010_04.CSV', 'NIMR_FEB2010_15.CSV', 'NIMR_FEB2010_16.CSV', 'NIMR_Feb2010_16.csv',
                                      'NIMR-report-Feb2010_03.csv', 'NIMR-report-Feb2010_06.csv', 'NIMR-report-Feb2010_05.csv', 'NIMR-report-Sep2009_03.csv', 'NIMR-report-Feb2010_04.csv']
        self.new_different_date_format = set()
        self.fix = set()

    def upload(self, ftype='flat', preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to TDB")
        measurements = self.parse(ftype, **kwargs)
        print('Formatting documents for upload')
        self.format_measurements(measurements, **kwargs)
        measurements = self.filter(measurements)
        measurements = self.create_index(measurements)
        #self.adjust_tdb_strain_names_from_vdb(measurements)
        print('Total number of indexes', len(self.indexes), 'Total number of measurements', len(measurements))
        if not preview:
            self.upload_documents(self.table, measurements, index='index', **kwargs)
        else:
            print("Titer Measurements:")
            print(json.dumps(measurements[0], indent=1))
            print("Remove \"--preview\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def format_measurements(self, measurements, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        self.HI_ref_name_abbrev =self.define_strain_fixes(self.HI_ref_name_abbrev_fname)
        self.define_location_fixes("source-data/flu_fix_location_label.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        for meas in measurements:
            meas['virus_strain'], meas['original_virus_strain'] = self.fix_name(self.HI_fix_name(meas['virus_strain'], serum=False))
            meas['serum_strain'], meas['original_serum_strain'] = self.fix_name(self.HI_fix_name(meas['serum_strain'], serum=True))
            self.test_location(meas['virus_strain'])
            self.test_location(meas['serum_strain'])
            self.add_attributes(meas, **kwargs)
            self.format_subtype(meas)
            self.format_assay_type(meas)
            self.format_date(meas)
            self.format_passage(meas, 'serum_passage', 'serum_passage_category')
            self.format_passage(meas, 'virus_passage', 'virus_passage_category')
            self.format_ref(meas)
            self.format_serum_sample(meas)
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
        return measurements

    def adjust_tdb_strain_names_from_vdb(self, measurements, database='vdb'):
        '''
        Compare measurement strain names to vdb strain names to ensure downstream matching between measurements and sequences
        '''
        table = self.virus + "_viruses"
        print("Using " + database + "." + table + " to adjust strain names to sequence strains")
        vdb_strains = self.relaxed_keys(set(list(r.db(database).table(table).get_field('strain').run())), self.relax_name)
        matched_virus_strains = 0
        matched_serum_strains = 0
        for meas in measurements:
            meas['virus_strain'], matched_virus_strains = self.adjust_strain_name(meas['virus_strain'], vdb_strains, matched_virus_strains)
            meas['serum_strain'], matched_serum_strains = self.adjust_strain_name(meas['serum_strain'], vdb_strains, matched_serum_strains)
        print(str(matched_virus_strains) + " out of " + str(len(measurements)) + " virus strains were matched to a sequence in " + database + "." + table)
        print(str(matched_serum_strains) + " out of " + str(len(measurements)) + " serum strains were matched to a sequence in " + database + "." + table)

    def adjust_strain_name(self, name, vdb_strains, number_matched):
        '''
        Change strain name if matched to a relaxed strain name in vdb strains
        Return adjusted strain name and total count of matched measurements
        '''
        relax_strain = self.relax_name(name)
        if relax_strain in vdb_strains:
            number_matched += 1
            if vdb_strains[relax_strain] != name:
                return vdb_strains[relax_strain], number_matched
        return name, number_matched

    def HI_fix_name(self, name, serum):
        '''
        Canonicalize strain names to match with vdb
        '''
        # replace all accents with ? mark
        original_name = name.encode('ascii', 'replace')
        # Replace whole strain names
        name = self.replace_strain_name(original_name, self.fix_whole_name)
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

    def replace_strain_name(self, original_name, fixes={}):
        '''
        return the new strain name that will replace the original
        '''
        if original_name in fixes:
            return fixes[original_name]
        else:
            return original_name

    def test_location(self, strain):
        '''
        Determine that strains come from known locations, if not, print suggestion to add location to
        flu_fix_location_label.tsv.
        '''
        if isinstance(strain, basestring) and "/" in strain:
            location = strain.split('/')[1]
            if self.determine_location(location) is None:
                print("Couldn't determine location for this strain, consider adding to flu_fix_location_label.tsv", location, strain)

    def add_attributes(self, meas, host, **kwargs):
        '''
        Add attributes to titer measurements
        '''
        meas['virus'] = self.virus.lower()
        if hasattr(self, 'subtype'):
            meas['subtype'] = self.subtype.lower()
        meas['host'] = host.lower()
        meas['timestamp'] = self.rethink_io.get_upload_timestamp()
        meas['inclusion_date'] = self.rethink_io.get_upload_date()

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

        meas['assay_date'] = meas['date']

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

    def format_subtype(self, meas):
        '''
        Format subtype attribute
        '''
        if 'subtype' in meas:
            if meas['subtype'] == 'H3':
                meas['subtype'] = 'h3n2'
            elif meas['subtype'] == 'H1 swl':
                meas['subtype'] = 'h1n1pdm'
            elif meas['subtype'] == 'B vic':
                meas['subtype'] = 'vic'
            elif meas['subtype'] == 'B yam':
                meas['subtype'] = 'yam'
        else:
            meas['subtype'] = None

    def format_assay_type(self, meas):
        '''
        Format assay_type attribute
        '''
        if 'assay_type' in meas:
            if meas['assay_type'] == 'HI':
                meas['assay_type'] = 'hi'
            elif meas['assay_type'] == 'FRA':
                meas['assay_type'] = 'fra'
        else:
            meas['assay_type'] = None

    def format_serum_sample(self,meas):
        '''
        Format serum sample attribute for the measurements, so that there is
        not a specific organism used for HI assays.
        '''
        if 'ferret_id' in meas.keys():
            new_id = meas['ferret_id']
            meas['serum_id'] = new_id
            meas['serum_host'] = 'ferret'
            meas.pop('ferret_id',None)
        if 'serum_id' in meas.keys() and meas['serum_id'].lower().startswith('f'):
            meas['serum_host'] = 'ferret'
        elif 'serum_id' in meas.keys() and meas['serum_id'].lower().startswith('h'):
            meas['serum_host'] = 'human'


    def create_index(self,  measurements, output=False):
        '''
        create unique key for storage in rethinkdb
        do this using an md5 hash on index_fields
        '''
        for meas in measurements:
            index = [meas[field] for field in self.index_fields]
            meas['index'] = hashlib.md5("".join(index)).hexdigest()
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
#       print("Filtering out measurements whose serum strain is not paired with a ref or test virus strain ensuring proper formatting")
#       measurements = filter(lambda meas: meas['serum_strain'] in self.ref_virus_strains or meas['serum_strain'] in self.test_virus_strains, measurements)
        print("Filtering out measurements missing required fields")
        measurements = filter(lambda meas: self.rethink_io.check_required_attributes(meas, self.upload_fields, self.index_fields, output=True), measurements)
        print("Filtering out measurements with virus, serum strain names or titers not formatted correctly")
        measurements = filter(lambda meas: self.correct_strain_format(meas['virus_strain'], meas['original_virus_strain']) and self.correct_strain_format(meas['serum_strain'], meas['original_serum_strain']) and self.correct_titer_format(meas['titer']), measurements)
        print(len(measurements), " measurements after filtering")
        return measurements

    def correct_titer_format(self, titer):
        if re.match("^[0-9<>.]*$", titer):
            return True
        else:
            print("Bad titer: {}. Removing.".format(titer))
            return False

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
