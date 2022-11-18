import os, re, time, datetime, csv, sys, json
from upload import upload
from rethinkdb import r
from Bio import SeqIO
import argparse
from parse import parse
from upload import parser
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload

parser.add_argument("--rename", action='store_true')

class cdc_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.removal_fields  = ['tested_by_fra', 'reported_by_fra', 'date', 'virus_collection_date', 'ref', 'virus_harvest_date', 'Boosted', 'RBC_species']
        self.cleanup_fields =  {'assay-type': 'assay_type', 'lot #': 'lot_number'}

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
        self.HI_ref_name_abbrev =self.define_strain_fixes(self.HI_ref_name_abbrev_fname)
        self.define_location_label_fixes("source-data/flu_fix_location_label.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        for meas in measurements:
            meas['source'] = 'cdc'
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
            self.format_titer(meas)
            if meas['ref'] == True:
                self.ref_serum_strains.add(meas['serum_strain'])
                self.ref_virus_strains.add(meas['virus_strain'])
            if meas['ref'] == False:
                self.test_virus_strains.add(meas['virus_strain'])
            meas['serum_host'] = self.get_serum_host(meas)
            # Update serum passage category for human pool serum using the serum id
            if meas['serum_host'] == 'human':
                self.format_passage(meas, 'serum_id', 'serum_passage_category')
            self.rethink_io.check_optional_attributes(meas, self.optional_fields)
            self.remove_fields(meas)
        if len(self.new_different_date_format) > 0:
            print("Found files that had a different date format, need to add to self.different_date_format")
            print(self.new_different_date_format)
        self.check_strain_names(measurements)
        return measurements

    def format_titer(self, meas):
        '''
        Format titer attribute.
        Correct "5" to "<20" in CDC uploads.
        '''
        if 'titer' in meas:
            if 'source' in meas:
                if meas['titer'] == '5' and meas['source'] == 'CDC':
                    meas['titer'] = '<20'
        else:
            meas['titer'] = None

    def format_serum_sample(self,meas):
        '''
        Format serum sample attribute for the measurements, so that there is
        not a specific organism used for HI assays.
        '''
        if not meas['serum_id']:
            if 'lot_number' in meas.keys():
                meas['serum_id'] = 'L' + str(meas['lot_number'])
            #     print('new serum id is %s' % (meas['serum_id']))
            # else:
            #     print('old serum id is %s' % (meas['serum_id']))

    def remove_fields(self, meas):
        '''
        Remove unnecessary fields provided in CDC titer tables.
        Specify fields for removal in self.removal_fields.
        '''
        for f in self.removal_fields:
            if f in meas.keys():
                meas.pop(f,None)

    def get_serum_host(self, meas):
        """
        Returns the serum host based on whether the host is in the measurement's serum id.
        The default serum host is 'ferret' if there are no matches for other hosts.
        """
        # These are known serum hosts that exist in the CDC titer data
        known_serum_hosts = ['human', 'mouse']
        for host in known_serum_hosts:
            if re.search(host, meas['serum_id'], re.IGNORECASE):
                return host

        return 'ferret'


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


def rename_cdc_columns(path, fstem) -> str:
    """
    Rename the columns according to previous CDC titer tables by using the
    column map file at `source-data/cdc_titer_column_map.tsv`. Assumes the
    titer file ends with .tsv and prints the filtered titers to a new TSV file.

    Returns the new reportable titer fstem which is the provided
    *fstem* + '_renamed'.
    """
    import pandas as pd

    column_map = pd.read_csv('source-data/cdc_titer_column_map.tsv', sep='\t', index_col="label")["fix"].to_dict()

    all_titer_file = path + fstem + '.tsv'
    all_titers = pd.read_csv(all_titer_file, sep='\t')\
                   .rename(columns=column_map)

    all_titers['assay_type'] = all_titers['assay_type'].str.replace('_protocol', '')

    renamed_titer_fstem = fstem + '_renamed'
    renamed_titer_file = path + renamed_titer_fstem + '.tsv'
    all_titers[column_map.values()].to_csv(renamed_titer_file, sep='\t', index=False)

    return renamed_titer_fstem


if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    if args.rename:
        args.fstem = rename_cdc_columns(args.path, args.fstem)

    connTDB = cdc_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
