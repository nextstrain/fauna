import os, re, time, datetime, csv, sys, json
from upload import upload
import rethinkdb as r
from Bio import SeqIO
import argparse
from parse import parse
from upload import parser
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload

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
        measurements = self.clean_field_names(measurements)        
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
            self.format_subtype(meas)
            self.format_assay_type(meas)                     
            meas['date'] = meas['assay_date']
            self.format_date(meas)
            self.format_passage(meas, 'serum_antigen_passage', 'serum_passage_category')
            self.format_passage(meas, 'virus_strain_passage', 'virus_passage_category')
            meas['serum_passage'] = meas.pop('serum_antigen_passage')
            meas['virus_passage'] = meas.pop('virus_strain_passage')
            meas.pop('passage',None)
            self.format_ref(meas)
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

    def remove_fields(self, meas):
        '''
        Remove unnecessary fields provided in CDC titer tables.
        Specify fields for removal in self.removal_fields.
        '''
        for f in self.removal_fields:
            if f in meas.keys():
                meas.pop(f,None)

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
