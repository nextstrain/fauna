import os, re, time, datetime, csv, sys, json
from upload import upload
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
parser.add_argument('--fstem', help="input file stem")
parser.add_argument('--ftype', default='flat', help="input file format, default \"flat\", other is \"tables\"")
parser.add_argument('--overwrite', default=False, action="store_true",  help ="Overwrite fields that are not none")
parser.add_argument('--exclusive', default=True, action="store_false",  help ="download all docs in db to check before upload")
parser.add_argument('--replace', default=False, action="store_true",  help ="If included, delete all documents in table")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('--preview', default=False, action="store_true",  help ="If included, preview a virus document to be uploaded")

class cdc_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.removal_fields  = ['tested_by_fra', 'reported_by_fra', 'date', 'virus_collection_date', 'ref']
        self.cleanup_fields =  {'assay-type': 'assay_type'}

    def upload(self, ftype='flat', preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to TDB")
        measurements = self.parse(ftype, **kwargs)
        print('Formatting documents for upload')
        self.format_measurements(measurements, **kwargs)
        measurements = self.clean_field_names(measurements)
        measurements = self.filter(measurements)
        measurements = self.create_index(measurements)
        self.adjust_tdb_strain_names(measurements)
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
        self.fix_whole_name.update(self.define_strain_fixes(self.HI_strain_fix_fname))
        self.HI_ref_name_abbrev =self.define_strain_fixes(self.HI_ref_name_abbrev_fname)
        self.define_location_fixes("source-data/flu_fix_location_label.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        for meas in measurements:
            meas['source'] = 'CDC'
            meas['virus_strain'], meas['original_virus_strain'] = self.fix_name(self.HI_fix_name(meas['virus_strain'], serum=False))
            meas['serum_strain'], meas['original_serum_strain'] = self.fix_name(self.HI_fix_name(meas['serum_strain'], serum=True))
            self.test_location(meas['virus_strain'])
            self.test_location(meas['serum_strain'])
            self.add_attributes(meas, **kwargs)
            self.format_date(meas)
            self.format_passage(meas, 'serum_antigen_passage', 'serum_passage_category')
            self.format_passage(meas, 'virus_strain_passage', 'virus_passage_category')
            meas['serum_passage'] = meas.pop('serum_antigen_passage')
            meas['virus_passage'] = meas.pop('virus_strain_passage')
            meas.pop('passage',None)
            self.format_id(meas)
            self.format_ref(meas)
            self.format_titer(meas)
            self.format_serum_sample(meas)
            if meas['ref'] == True:
                self.ref_serum_strains.add(meas['serum_strain'])
                self.ref_virus_strains.add(meas['virus_strain'])
            if meas['ref'] == False:
                self.test_virus_strains.add(meas['virus_strain'])
            self.rethink_io.check_optional_attributes(meas, self.optional_fields)
            self.remove_fields(meas)
        if len(self.new_different_date_format) > 0:
            print("Found files that had a different date format, need to add to self.different_date_format")
            print(self.new_different_date_format)
        self.check_strain_names(measurements)
        return measurements

    def format_date(self, meas):
        '''
        Format date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        Input date could be YYYY_MM_DD, reformat to YYYY-MM-DD
        '''

        meas['date'] = meas['assay_date']
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

    def remove_fields(self, meas):
        '''
        Remove unnecessary fields provided in CDC titer tables.
        Specify fields for removal in self.removal_fields.
        '''
        for f in self.removal_fields:
            if f in meas.keys():
                meas.pop(f,None)
                print f

    def clean_field_names(self, measurements):
        '''
        Change field names from CDC titer tables to fit into fauna db.
        Dictionary of field names to update {'old_name': 'new_name'} stored in self.cleanup_fields.
        '''
        for meas in measurements:
            for old_name in self.cleanup_fields.keys():
                if old_name in meas.keys():
                    new_name = self.cleanup_fields[old_name]
                    meas[new_name] = meas[old_name]
                    meas.pop(old_name,None)
        return measurements

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = cdc_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
