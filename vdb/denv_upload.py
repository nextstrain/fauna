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


    def fix_name(self, doc): # this finally gets called in the last line, which calls upload (method in upload.py), which calls format_viruses.
        original_name = doc['strain'] # keep copy of original name

        if self.replace_strain_name(original_name, self.fix_whole_name) != original_name:
            name = (self.replace_strain_name(original_name, self.fix_whole_name)).lower() # if we've manually edited this specific name, leave that alone.
        else: # Otherwise, make a new standard name.
            s = None # First try and identify serotype from direct annotation or from the strain name.
            if len(doc['serotype'].split('_'))==3: # dengue_virus_1234 --> 1234
                s = doc['serotype'].split('_')[2]

            else: # some version of D1/DENV1/DEN1 annotation in name
                patterns = ['[\/\\_-]?DEN[1-4]{1}[\/\\_-]?', '[\/\\_-]?DENV[1-4]{1}[\/\\_-]?', '[\/\\_-]?D[1-4]{1}[\/\\_-]?']
                for p in patterns: # Try each pattern, break if we find one that matches.
                    try:
                        s = re.search(p, original_name).group(0) # /DENV1/
                        s = re.search('[1-4]{1}', s).group(0) # isolate the serotype number
                        break
                    except:
                        continue
            try:
                sero = 'd'+str(int(s)) # check if we isolated an integer, return as string
            except: #
                # print "Strain %s has apparent serotype %s, cannot parse; setting to `d`"%(original_name, doc['serotype'])
                sero = 'd'

            try: # try to parse the year from the collection date
                yr = str(int(doc['collection_date'][:4])) # Year is usually first 4 digits of collection date
                assert yr[:2] == '19' or yr[:2]=='20'
                name = '/'.join([sero, doc['country'], yr]).lower() # New name: d1234/country/year
            except:
                # print 'No year found for strain %s, setting to None'%original_name
                # print 'Collection date annotation:', doc['collection_date']
                yr = None
                name = ('%s/%s/None'%(sero, doc['country'])).lower()
        return name, original_name


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
