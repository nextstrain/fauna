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

parser.add_argument('--assay_type', default='HI', help='type of assay being recorded')

class nimr_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.assay_type = kwargs['assay_type']
        self.assay_date = None

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
            meas['assay_type'] = self.assay_type
            self.test_location(meas['virus_strain'])
            self.test_location(meas['serum_strain'])
            self.add_attributes(meas, **kwargs)
            self.format_subtype(meas)
            self.format_assay_type(meas)                
            self.format_date(meas)
            if meas['assay_date'] != None: self.assay_date = meas['assay_date']
            self.format_passages(meas)
            self.format_ref(meas)
            self.format_serum_sample(meas)
            if meas['ref'] == True:
                self.ref_serum_strains.add(meas['serum_strain'])
                self.ref_virus_strains.add(meas['virus_strain'])
            if meas['ref'] == False:
                self.test_virus_strains.add(meas['virus_strain'])
            self.rethink_io.check_optional_attributes(meas, self.optional_fields)
            meas['assay_type'] = self.assay_type
        if len(self.new_different_date_format) > 0:
            print("Found files that had a different date format, need to add to self.different_date_format")
            print(self.new_different_date_format)
        self.check_strain_names(measurements)
        for meas in measurements:
            meas['assay_date'] = self.assay_date
        return measurements

    def format_passages(self, meas, source_type='NIMR'):
        '''
        Create virus passage and category fields for NIMR documents and
        null fields for serum passage and category
        '''
        self.format_passage(meas, 'passage', 'passage_category')
        if source_type == 'NIMR':
            meas['virus_passage'] = meas['passage']
            meas['virus_passage_category'] = meas['passage_category']
            meas['serum_passage'] = False
            meas['serum_passage_category'] = False
            meas.pop('passage',None)
            meas.pop('passage_category',None)
        else:
            meas['serum_passage'] = meas['passage']
            meas['serum_passage_category'] = meas['passage_category']
            meas['virus_passage'] = False
            meas['virus_passage_category'] = False
            meas.pop('passage',None)
            meas.pop('passage_category',None)

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = nimr_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
