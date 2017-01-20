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

class elife_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

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
        self.check_uniqueness(measurements)
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
            meas['virus_strain'], meas['original_virus_strain'] = self.fix_name(self.HI_fix_name(meas['virus_strain'], serum=False))
            meas['serum_strain'], meas['original_serum_strain'] = self.fix_name(self.HI_fix_name(meas['serum_strain'], serum=True))
            self.test_location(meas['virus_strain'])
            self.test_location(meas['serum_strain'])
            self.add_attributes(meas, **kwargs)
            self.format_date(meas)
            meas['assay_date'] = "XXXX-XX-XX"
            self.format_passage(meas, 'serum_passage', 'serum_passage_category')
            self.format_passage(meas, 'virus_passage', 'virus_passage_category')
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
        if len(self.new_different_date_format) > 0:
            print("Found files that had a different date format, need to add to self.different_date_format")
            print(self.new_different_date_format)
        self.check_strain_names(measurements)
        self.disambiguate_sources(measurements)
        return measurements

    def disambiguate_sources(self, measurements):
        sources = {}
        for meas in measurements:
            src = meas['source']
            if src not in sources.keys():
                sources[src] = 0
            else:
                sources[src] += 1
            new_src = src + '_' + str(sources[src])
            meas['source'] = new_src

    def check_uniqueness(self, measurements):
        indices = []
        unique = 0
        nonunique = 0
        uniq = True
        for meas in measurements:
            index_string = ''
            for field in meas['index']:
                index_string = index_string + str(field)
            if index_string in indices:
                print "Nonunique index field: ", index_string
                nonunique += 1
                uniq = False
            else:
                indices.append(index_string)
                unique += 1
        print "Unique fields: ", unique
        print "Nonunique fields: ", nonunique
        return uniq

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = elife_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
