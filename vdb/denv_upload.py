import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class denv_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs) # sets virus, viruses_table, sequences_table, database, uploadable_databases;
                                        # initiates empty objects for strains, strain_fix_fname, virus_to_sequence_transfer_fields
                                        # also calls parse.py __init__ --> set gbdb, checks for accession-type data in kwargs.

        self.strain_fix_fname = "source-data/denv_strain_name_fix.tsv" # tsv with label\tfix names

    def format_viruses(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        if self.strain_fix_fname is not None:
            self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        self.define_regions("source-data/geo_regions.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        self.define_latitude_longitude("source-data/geo_lat_long.tsv", "source-data/geo_ISO_code.tsv")
        for doc in documents:
            self.format_date(doc)
            self.format_place(doc)
            self.format_region(doc)
            self.determine_latitude_longitude(doc, ['country'])
            doc['strain'], doc['original_strain'] = self.fix_name(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)

    def format_sequences(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        if self.strain_fix_fname is not None:
            self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        for doc in documents:
            self.format_place(doc)
            self.format_date(doc)
            doc['strain'], doc['original_strain'] = self.fix_name(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)


    def fix_name(self, doc):
        '''
        Given metadata annotations in doc, make new strain names like
        DENV1234/country/ID/year
        where ID is derived from the original strain ID with metadata redundancies removed
        '''

        original_name = doc['strain'] # keep copy of original name

        try: # try to parse the year from the collection date
            yr = str(int(doc['collection_date'][:4])) # Year is usually first 4 digits of collection date
        except:
            yr = 'NA'

        try:
            s = 'DENV'+str(int(doc['serotype'].split('_')[2])) # Check for vipr serotype annotation like 'dengue_virus_1234'
        except:
            s = 'DENV'

        country, division, region = doc['country'], doc['division'], doc['region']
        try:
            country_code = self.country_to_code[country]
        except:
            country_code='NA'

        not_metadata = [] # Pull the parts of the original name that are not metadata redundancies.

        split_name = re.split('[\W^_]+', original_name) # Split on any punctuation, look at each element individually.

        for element in split_name: # Check each element to see if it looks like common kinds of metadata.
            is_metadata = False

            year_patterns = ['^%s$'%yr, '^%s$'%yr[-2:]]
            geo_patterns = [country, division, region, '^%s$'%country_code ] ## Do the same thing we did for year, but with the country, region, city, and country code.
            strain_patterns = ['DEN[1-4]{1}', 'DENV[1-4]{1}', 'D[1-4]{1}']

            for p in year_patterns+geo_patterns: # Once we identify an element as metadata, we can chuck it.
                if re.search(p, element, re.IGNORECASE):
                    is_metadata = True
                    break

            if s == 'DENV' or is_metadata==False: # We still want to look for serotype annotations that got missed by vipr.
                for p in strain_patterns:
                    try:
                        s = re.search(p, element).group(0) # DENV1
                        s = 'DENV'+re.search('[1-4]{1}', s, re.IGNORECASE).group(0) # Reformat --> 'DENV1234'
                        is_metadata = True # flag the element as metadata
                        break
                    except:
                        continue
            if is_metadata == False:
                not_metadata.append(element)

        if not_metadata == []:
            new_id = 'NA' # figure out what to make the new id when they don't give one...
        else:
            new_id = '_'.join(not_metadata)

        name = '%s/%s/%s/%s'%(s, doc['country'], new_id, yr)

        print original_name, new_id
        return original_name, name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    args = parser.parse_args() # parser is an argparse object initiated in parse.py
    virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country', 7:'serotype'} # define fields in fasta headers that you want used in parse.py > parse > parse_fasta_file ---> (viruses, sequences)
    sequence_fasta_fields = {0:'accession', 1:'strain', 2: 'locus', 3: 'collection_date', 5: 'country', 7: 'serotype'}
    # 0         1   2   3       4       5   6   7
    #>KP780172|#39|NA|2014_09_30|Human|China|NA|Dengue_virus_1234

    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = denv_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
