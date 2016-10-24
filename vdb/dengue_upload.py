import os, re, time, datetime, csv, sys, json
import rethinkdb as r
import pandas as pd
from Bio import SeqIO
from upload import upload
from upload import parser

class dengue_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs) # sets virus, viruses_table, sequences_table, database, uploadable_databases;
                                        # initiates empty objects for strains, strain_fix_fname, virus_to_sequence_transfer_fields
                                        # also calls parse.py __init__ --> set gbdb, checks for accession-type data in kwargs.
        self.strain_fix_fname = "source-data/dengue_strain_name_fix.tsv" # tsv with label\tfix names

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        self.connect(**kwargs)
        print("Uploading Viruses to VDB")
        data = self.parse(**kwargs)[0]
        print('Formatting documents for upload')
        self.format_metadata(data, **kwargs)
        viruses, sequences = self.separate_viruses_sequences(data, **kwargs)
        print("")
        print("Filtering out viruses")
        viruses = self.filter(viruses, 'strain', **kwargs)
        print("Filtering out sequences")
        sequences = self.filter(sequences, 'accession', **kwargs)
        print("")
        self.match_duplicate_strains(viruses, sequences, **kwargs)
        self.match_database_duplicate_strains(viruses, sequences, **kwargs)
        self.match_duplicate_accessions(sequences, **kwargs)
        self.match_database_duplicate_accessions(sequences, **kwargs)
        self.link_viruses_to_sequences(viruses, sequences)
        #self.transfer_fields(viruses, sequences, self.virus_to_sequence_transfer_fields)
        print("")
        print("Upload Step")
        if not preview:
            print("Uploading viruses to " + self.database + "." + self.viruses_table)
            self.upload_documents(self.viruses_table, viruses, index='strain', **kwargs)
            print("Uploading sequences to " + self.database + "." + self.sequences_table)
            self.upload_documents(self.sequences_table, sequences, index='accession', **kwargs)
        else:
            print("Viruses:")
            print(json.dumps(viruses[0], indent=1))
            print("Sequences:")
            print(json.dumps(sequences[0], indent=1))
            print("Remove \"--preview\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def parse_tsv_file(self, tsv, **kwargs):
        data = pd.read_csv(tsv, dtype='str', delimiter='\t', index_col=0, header=0, skiprows=1)# Read in tsv --> dataframe
        header_rename = {'Accession': 'accession', 'Name': 'strain', 'Segment': 'locus', 'Start': 'start', 'Stop': 'stop', 'Country': 'country', 'Sampling Year': 'collection_date', 'Organism': 'serotype', 'Species': 'virus', 'Isolate Name': 'isolate_name', 'Georegion': 'region', 'Author': 'authors', 'Sampling City': 'location', 'Pubmed ID': 'PMID', 'Sequence': 'sequence' }
        data.rename(columns=header_rename, inplace=True) # rename columns
        data.drop([i for i in ['code','Patient Id', 'Sequence Length','ssam_tx_gb_taxid','ssam_tx_full_lineage','Genus','Family'] if i in data.columns.values],inplace=True,axis=1) # drop superfluous columns
        data = data.where((pd.notnull(data)), None) # Fill np.nan floats with nonetype objects to maintain compatibility
        data = data.to_dict(orient='index').values() # dataframe --> [ {column_name: value_row_0, ...}, ... {column_name: value_row_nsequences, ...} ]
        return data, ''

    def format_metadata(self, documents, **kwargs):
        if self.strain_fix_fname is not None:
            self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        self.define_regions("source-data/geo_regions.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        self.define_latitude_longitude("source-data/geo_lat_long.tsv", "source-data/geo_ISO_code.tsv")
        for doc in documents:
            self.format_date(doc)
            if doc['location'] != None:
                doc['location'] = re.sub(r'[\s\W]','', doc['location'])
                if doc['location'].lower() in ['south', 'northern', 'centralwest', 'westernprovince', 'northeast']:
                    doc['location'] = None
            if doc['country'] != None:
                doc['country'] = re.sub(r'[\s\W]','', doc['country'])
                if doc['country'].lower() in ['south', 'northern', 'centralwest', 'westernprovince']:
                    doc['country'] = None
            self.format_place(doc)
            self.format_region(doc)
            self.determine_latitude_longitude(doc, ['location', 'country'])
            self.fix_strain(doc)
            self.fix_locus(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)

    def fix_locus(self, doc, **kwargs):
        reference_loci = {'C': (95, 436), 'prM': (437, 934), 'E': (935, 2413), 'NS1': (2414, 3469), 'NS2A': (3470, 4123), 'NS2B': (4124, 4513), 'NS3': (4514, 6370), 'NS4A': (6371, 6820), 'NS4B': (6821, 7564), 'NS5': (7565, 10264)}
        seq_loci = []
        seq_start, seq_end = int(doc['start']), int(doc['stop'])
        for locus, (ref_start, ref_end) in reference_loci.items():
            match_start = max(ref_start, seq_start)
            match_end = min(ref_end, seq_end)
            locus_length = float(ref_end-ref_start)
            prop_match = float(match_end - match_start) / locus_length # < 0 if entirely outside of the reference coordinates.
            if prop_match >= 0.8: # call it a match if >=80% of the locus is covered by the sequence.
                seq_loci.append(locus)
        doc['locus'] = seq_loci

    def separate_viruses_sequences(self, data, **kwargs):
        viruses = []
        sequences = []
        for record in data:
            v = {k: v for k,v in record.items() if k in virus_attribs}
            s = {k: v for k,v in record.items() if k in sequence_attribs}
            v = self.add_virus_fields(v, **kwargs)
            s = self.add_sequence_fields(s, **kwargs)
            sequences.append(s)
            viruses.append(v)
        return (viruses, sequences)


    def fix_strain(self, doc, **kwargs):
        '''
        Given metadata annotations in doc, make new strain names like
        DENV1234/country/ID/year
        where ID is derived from the original strain ID with metadata redundancies removed
        '''
        original_strain = doc['strain'] # Keep copy of original name

        if self.replace_strain_name(original_strain, self.fix_whole_name) != original_strain: # If we've manually edited this specific name, leave that alone.
            doc['strain'], doc['original_strain'] = (self.replace_strain_name(original_strain, self.fix_whole_name)).upper(), original_strain
            return
        elif doc['isolate_name']!=None and doc['isolate_name']!=original_strain: # Check for alternate ID
            strain = doc['isolate_name']
            print 'Accession %s has isolate name %s and strain name %s. Using isolate name.'%(doc['accession'], doc['isolate_name'], original_strain)
        else:
            strain = original_strain
        try:
            del doc['isolate_name']
        except KeyError:
            pass

        ### Pull metadata from documentation
        if doc['country'] != None:
            country = doc['country']
        else:
            country = 'NA'
        if doc['location'] != None:
            location = doc['location']
        else:
            location = 'NA'
        try:
            year = doc['collection_date'].split('-')[0]
        except AttributeError:
            year = 'NA'
        try:
            sero = 'DENV'+str(int(doc['serotype'][-1])) # Check for serotype annotation like 'dengue virus 1234'
        except ValueError:
            sero = 'NA'
        try:
            country_code = self.country_to_code[country]
        except KeyError:
            country_code = 'NA'

        ### Remove metadata redundancies from the strain ID
        # date_patterns = [r'%s[\W^_]?%s[\W^_]?%s'%(day, month, year), r'%s[\W^_]?%s[\W^_]?%s'%(day, month, year[-2:]), # Full dates
        # r'%s[\W^_]?%s[\W^_]?%s'%(month, day, year), r'%s[\W^_]?%s[\W^_]?%s'%(month, day, year[-2:]),
        # r'%s[\W^_]?%s'%(month, year), r'%s[\W^_]?%s'%(month, year[-2:]), # month-year
        date_patterns = [r'%s'%year, r'Y%s(\b|_)'%year[-2:], r'(\b|_)Y%s'%year[-2:], r'(\b|_)%s'%year[-2:], r'%s(\b|_)'%year[-2:], #r'(\b|_)%s(\b|_)'%month,  Just month or just year (offset by punctuation)
        r'(\b|_)%s%s(\b|_)'%(year[-2:], country_code), r'(\b|_)%s%s(\b|_)'%(country_code, year[-2:]) # country code + year
        ]
        disease_patterns = ['DF', 'DHF', 'AS', 'DSS', 'U', 'UD']
        geo_patterns = [r'(\b|_)%s(\b|_)'%country_code, r'(\b|_)%s[\w]{1}(\b|_)'%country_code,
        r'(\b|_)[\w]{1}%s(\b|_)'%country_code, r'(\b|_)%s[\w]{1}%s(\b|_)'%(country_code[0], country_code[1]),] # 2 or 3-letter country code (alone)
        host_patterns = [r'human', r'unknown', r'mosquito']
        strain_patterns = [r'(\b|_)DF(\b|_)', r'DEN[\W^_]?[1-4]{1}', r'DENV[\W^_]?[1-4]{1}', r'D[\W^_]?[1-4]{1}'] # Check for serotype annotation in the metadata.
        all_patterns = disease_patterns+date_patterns+strain_patterns+geo_patterns+host_patterns

        for p in all_patterns:
            strain = re.sub(p, '', strain, re.IGNORECASE)
        strain = re.sub(r'[\W^_]', '', strain)                     # Remove punctuation and whitespace
        strain = strain.upper().replace(country.upper(), '')       # Last sweep for multi-word countries (e.g. East_Timor)
        strain = strain.upper().replace(location.upper(), '')       # Last sweep for multi-word cities

        if strain == '':
            strain = ''
        doc['original_strain'] = original_strain
        doc['strain'] = ('%s/%s/%s/%s'%(sero, country, strain, year)).strip().upper() # New strain name = DENV1234/COUNTRY/STRAIN_ID/YEAR
        doc['serotype'] = sero                                             # Update serotype in standard format

        print original_strain
        print doc['strain'], '\n\n'


    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    args = parser.parse_args() # parser is an argparse object initiated in parse.py
    virus_attribs = ['strain', 'original_strain', 'virus', 'serotype','collection_date', 'region', 'country', 'location','authors', 'PMID'] # define fields in fasta headers that you want used in parse.py > parse > parse_fasta_file ---> (viruses, sequences)
    sequence_attribs = ['accession', 'strain', 'original_strain', 'virus', 'serotype',  'locus', 'collection_date', 'sequence']

    setattr(args, 'virus_attribs', virus_attribs)
    setattr(args, 'sequence_attribs', sequence_attribs)
    connVDB = dengue_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
