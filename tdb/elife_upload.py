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
        fstem_assay_date = self.parse_assay_date_from_filename(kwargs['fstem'])
        for meas in measurements:
            meas['virus_strain'], meas['original_virus_strain'] = self.fix_name(self.HI_fix_name(meas['virus_strain'], serum=False))
            meas['serum_strain'], meas['original_serum_strain'] = self.fix_name(self.HI_fix_name(meas['serum_strain'], serum=True))
            self.test_location(meas['virus_strain'])
            self.test_location(meas['serum_strain'])
            self.add_attributes(meas, **kwargs)
            self.format_subtype(meas)
            self.format_assay_type(meas)
            self.format_date(meas)
            if meas.get('assay_date') is None:
                meas['assay_date'] = fstem_assay_date
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
        self.disambiguate_sources(measurements)
        return measurements


    def parse_assay_date_from_filename(self, fstem):
        """
        Parse assay date from the *fstem*.
        *fstem* is expected to be formatted as `YYYYMMDD*`

        If unable to parse date from *fstem*, then return masked assay date as `XXXX-XX-XX`.
        If there are multiple valid dates, then return the latest date.
        """
        assay_date = "XXXX-XX-XX"
        valid_dates = set()
        for potential_date in re.findall(r"\d{8}", fstem):
            # Check if the dates are valid
            try:
                date = datetime.datetime.strptime(potential_date, '%Y%m%d')
            except ValueError:
                continue

            # Date is only a valid assay date if it's earlier than the current datetime!
            if date < datetime.datetime.now():
                valid_dates.add(date)

        if len(valid_dates) == 0:
            print(f"Failed to parse assay date from filename {fstem!r}")
        elif len(valid_dates) == 1:
            assay_date = datetime.datetime.strftime(valid_dates.pop(), '%Y-%m-%d')
        else:
            sorted_dates = [datetime.datetime.strftime(valid_date, '%Y-%m-%d') for valid_date in sorted(valid_dates)]
            raise Exception(f"Found multiple potential assay dates in filename {fstem!r}: {sorted_dates}. " +
                            "Filename should only contain one assay date in format YYYYMMDD.")

        return assay_date


    def disambiguate_sources(self, measurements):
        '''
        Add counter to sources so that create_index still creates unique identifiers for each
        titer value.
        '''
        sources = {}
        for meas in measurements:
            src = meas['source']
            if src is None:
                src = 'UnknownSource'
            if src not in sources.keys():
                sources[src] = 0
            else:
                sources[src] += 1
            new_src = src + '_' + str(sources[src])
            meas['source'] = new_src

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = elife_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
