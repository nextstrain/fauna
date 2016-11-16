import os, re, time, datetime, csv, sys, json
import rethinkdb as r
import pandas as pd
from Bio import SeqIO
from upload import upload
from upload import parser
from update import update

class dengue_upload(update):
    def __init__(self, **kwargs):
        # upload.__init__(self, **kwargs) # update: sets virus, viruses_table, sequences_table, database, uploadable_databases;
                                        # initiates empty objects for strains, strain_fix_fname, virus_to_sequence_transfer_fields
                                        # parse: set gbdb, checks for accession-type data in kwargs.
        update.__init__(self, **kwargs) # update(upload(parse))
        self.nonames=[]

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        self.connect(**kwargs)
        print("Uploading Viruses to VDB")
        data = self.parse(**kwargs)[0]                      # Calls parse_tsv_file below
        self.remove_duplicate_sequences(data)               # Keeps oldest record of identical accession/sequence pairs. Also sorts by accession number.
        print('Formatting documents for upload')
        self.format_metadata(data, **kwargs)                # Pull source data, format place, date, strain, and citation metadata
        viruses, sequences = self.separate_viruses_sequences(data, **kwargs) # Split into separate virus and sequence objects,
                                                                             # with multiple records from the same isolate as separate
                                                                             # sequences linked to the same virus object
        print("\nFiltering out viruses")
        viruses = self.filter(viruses, 'strain', **kwargs)
        print("Filtering out sequences\n")
        sequences = self.filter(sequences, 'accession', **kwargs)

        self.match_duplicate_strains(viruses, sequences, **kwargs)
        self.match_database_duplicate_strains(viruses, sequences, **kwargs)
        self.match_duplicate_accessions(sequences, **kwargs)
        self.match_database_duplicate_accessions(sequences, **kwargs)
        self.link_viruses_to_sequences(viruses, sequences)

        print 'No name strains:'
        print self.nonames

        print("\nUpload Step")
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
        return data, '' # extra '' purely for syntactical reasons; ignored later.

    def format_metadata(self, documents, **kwargs):
        self.define_regions("source-data/geo_regions.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        self.define_latitude_longitude("source-data/geo_lat_long.tsv", "source-data/geo_ISO_code.tsv")
        self.get_genbank_dates(documents, **kwargs)
        for doc in documents:
            try:
                doc['serotype'] = str(int(doc['serotype'][-1]))
            except:
                doc['serotype'] = None
            self.format_date(doc) # overriden below
            self.format_place(doc) # overriden below
            self.format_region(doc)
            self.determine_latitude_longitude(doc, ['location', 'country'])
            self.fix_strain(doc)   # overridden below
            try:
                doc['authors'] = doc['authors'].split(',')[0].split(' ')[-1]+' et al.' # just keep the first author
            except:
                pass
            self.fix_locus_gene_list(doc)
            self.rethink_io.check_optional_attributes(doc, [])

    def get_genbank_dates(self, documents, **kwargs):
        '''
        Update all fields using genbank files
        '''
        print("Updating collection dates from genbank")
        accessions = ','.join(list(set([ k['accession'].strip() for k in documents ])))
        self.accessions=accessions
        dates = self.get_genbank_sequences(**kwargs)[0]
        dates = { k:v for (k,v) in dates }
        for doc in documents:
            if doc['accession'] in dates:
                doc['collection_date'] = dates[doc['accession']]
        print 'Updated genbank dates for %d documents'%len(dates)

    def parse_gb_entries(self, handle, **kwargs):
        '''
        Go through genbank records to get proper collection dates
        '''
        dates = []
        SeqIO_records = SeqIO.parse(handle, "genbank")
        for record in SeqIO_records:
            accession = re.match(r'^([^.]*)', record.id).group(0).upper()  # get everything before the '.'?
            for feat in record.features:
                if feat.type == 'source' and 'collection_date' in feat.qualifiers:
                    dates.append( (accession, self.convert_gb_date(feat.qualifiers['collection_date'][0]) ) )
        handle.close()
        print("There were " + str(len(dates)) + " records in the parsed file")
        return (dates, [])

    def convert_gb_date(self, collection_date):
        '''
        Converts calendar dates between given formats
        '''
        N_fields = len(collection_date.split('-'))
        if N_fields == 1:
            return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%Y'), '%Y-XX-XX')
        elif N_fields == 2:
            try:
                return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%b-%Y'), '%Y-%m-XX')
            except:
                return '%s-XX'%(collection_date)
        elif N_fields == 3:
            try:
                return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%d-%b-%Y'), '%Y-%m-%d')
            except:
                return collection_date

    def remove_duplicate_sequences(self, documents, **kwargs):
        accessions = list(set([ doc['accession'] for doc in documents ])) # all accessions in data
        accessions_to_documents = { a:[doc for doc in documents if doc['accession'] == a] for a in accessions}
        # { accession:[record_dictionary1, record_dictionary2, ...]}
        number_duplicates = 0
        for a, docs in accessions_to_documents.iteritems():
            sequences = set([ doc['sequence'].lower() for doc in docs ]) # check all the sequences for this accession are identical
            if len(sequences) > 0 or len(docs)<2:
                continue

            docs.sort(key=lambda doc: doc['collection_date']) # sort records by collection date
            oldest = docs[0]                                  # keep the oldest, delete the rest
            duplicates = docs[1:]
            for d in duplicates:
                documents.remove(d)
                number_duplicates += 1
        documents.sort(key=lambda doc: doc['accession']) # sort by accession number to make strain name debugging easier later
        if number_duplicates > 0:
            print 'Removed %d duplicate accessions with the same sequence, kept the oldest record.'%number_duplicates

    def format_place(self, doc, determine_location=True):
        '''
        Try to determine location information from location fields
        Ensure snakecase formatting after assigning location fields
        '''
        location_fields = ['location', 'division', 'country']
        for field in location_fields:
            if determine_location: # Check against geo source data?
                if field in doc and doc[field] is not None:
                    # Cleanups specific to LANL formatting of dengue metadata
                    doc[field] = re.sub(r'[\s\W]','', doc[field])
                    doc[field] = doc[field].lower().replace('ningbodiseasepreventionandcontrolcenterzhejiangprovince', 'Ningbo')
                    if doc[field].lower().endswith('province'):
                        doc[field] = doc[field].lower().replace('province', '')
                    if field in ['country', 'location'] and doc[field].lower() in ['south', 'northern', 'centralwest', 'western', 'northeast']:
                        doc[field] = None
                        continue

                    # Check against synonymns in geo source data
                    result = self.determine_location(self.snakecase_to_camelcase(doc[field]))
                    if result is not None:
                        doc['location'], doc['division'], doc['country'] = result
                        break
                    else:
                        print("couldn't parse location for ", doc['strain'], self.snakecase_to_camelcase(doc[field]))
                        if "_" in doc[field]:  # French_Polynesia -> french_polynesia
                            doc[field] = "_".join(doc[field].split("_")).lower()
                        doc[field] = self.camelcase_to_snakecase(doc[field])
            else:
                if field in doc and doc[field] is not None:
                    if "_" in doc[field]:  # French_Polynesia -> french_polynesia
                        doc[field] = "_".join(doc[field].split("_")).lower()
                    doc[field] = self.camelcase_to_snakecase(doc[field])

    def fix_locus_gene_list(self, doc, **kwargs):
        reference_loci = {'C': (95, 436), 'prM': (437, 934), 'E': (935, 2413), 'NS1': (2414, 3469), 'NS2A': (3470, 4123), 'NS2B': (4124, 4513), 'NS3': (4514, 6370), 'NS4A': (6371, 6820), 'NS4B': (6821, 7564), 'NS5': (7565, 10264)}
        seq_loci = []
        seq_start, seq_end = int(doc['start']), int(doc['stop']) # coordinates of this sequence
        for locus, (ref_start, ref_end) in reference_loci.items(): # coordinates of the locus
            match_start = max(ref_start, seq_start)
            match_end = min(ref_end, seq_end)
            locus_length = float(ref_end-ref_start)
            prop_match = float(match_end - match_start) / locus_length # < 0 if entirely outside of the reference coordinates.
            if prop_match >= 0.8: # call it a match if >=80% of the locus is covered by the sequence.
                seq_loci.append(locus)
        if seq_loci > 1: # if captured more than one gene, locus just == 'genome', keep the full list in an optional field (gene_list)
            doc['locus'] = 'genome'
            doc['gene_list'] = seq_loci
        else:
            doc['locus'] = seq_loci[0] # if only one gene in the sequence, both locus and the gene_list are the gene
            doc['gene_list'] = seq_loci[0]

    def fix_strain(self, doc, **kwargs):
        '''
        Given metadata annotations in doc, make new strain names like
        DENV1234/country/ID/year
        where ID is derived from the original strain ID with metadata redundancies removed
        '''
        ## Pull data, try the obvious first
        strain_id = doc['strain']
        doc['original_strain'] = doc['strain'] # Keep copy of original name
        country_code, country, division, location, sero, year, month, day = self.pull_metadata(doc) # Pull metadata from pre-processed annotations

        if doc['isolate_name']!=None and doc['isolate_name']!=strain_id:    # Check for alternate ID field
            strain_id = doc['isolate_name']
            print 'Accession %s has isolate name %s and strain name %s. Using isolate name.'%(doc['accession'], doc['isolate_name'], original_strain)
        try:                                                             # Many sequences already have a Broad ID. If so, use it.
            strain_id = re.search(r'BID[-]?V?[0-9]+', strain_id).group(0)   # e.g. BID-V123456, - and/or V optional.
            strain_id = re.sub(r'[\W^_]', '', strain_id)
            doc['strain'] = self.build_canonical_name(sero, country, strain_id, year)
            return
        except:
            pass

        ### Remove metadata redundancies from the strain ID
        disease_patterns = ['DF', 'DHF', 'AS', 'DSS', 'UD', 'EHI',r'(\b|_)NA(\b|_)']
        date_patterns = [r'%s[\W^_]?%s[\W^_]?%s'%(day, month, year), r'%s[\W^_]?%s[\W^_]?%s'%(day, month, year[-2:]), # Full dates
        r'%s[\W^_]?%s[\W^_]?%s'%(month, day, year), r'%s[\W^_]?%s[\W^_]?%s'%(month, day, year[-2:]),
        r'%s[\W^_]?%s'%(month, year), r'%s[\W^_]?%s'%(month, year[-2:]),                                              # month-year
        r'%s'%year, r'M?Y%s(\b|_)'%year[-2:], r'(\b|_)M?Y%s'%year[-2:], r'(\b|_)%s'%year[-2:], r'%s(\b|_)'%year[-2:], # just year (offset by punctuation)
        r'(\b|_)%s%s(\b|_)'%(year[-2:], country_code), r'(\b|_)%s%s(\b|_)'%(country_code, year[-2:])] # country code + year
        geo_patterns = [r'(\b|_)%s'%country_code, r'%s(\b|_)'%country_code, r'(\b|_)%s[\w]{1}(\b|_)'%country_code, # 2 or 3-letter country code (alone)
        r'(\b|_)[\w]{1}%s(\b|_)'%country_code, r'(\b|_)%s[\w]{1}%s(\b|_)'%(country_code[0], country_code[1])]
        host_patterns = [r'human', r'unknown', r'mosquito', r'(\b|_)h(\b|_)', r'hu(\b|_)', r'(\b|_)hu', r'homosapien', r'homosapiens']
        strain_patterns = [r'DEN[\W^_]?[1-4]{1}', r'DENV[\W^_]?[1-4]{1}', r'(\b|_)DV?[\W^_]?[1-4]{1}'] # serotype annotations
        punctuation_and_whitespace = [r'[\W^_]'] # word boundaries are useful above, but we want to remove them eventually.
        multi_word_patterns = [country, location, division] # these are often multi-word; best to look after removing word boundaries.
        all_patterns = disease_patterns+date_patterns+strain_patterns+geo_patterns+host_patterns # Loop through twice because sometimes order matters.
        all_patterns = all_patterns*2 + punctuation_and_whitespace+multi_word_patterns           # Remove the punctuation and multi-word patterns once at the end.

        for p in all_patterns: # Try and remove each of the above patterns.
            strain_id = re.sub(p, '', strain_id.lower(), flags=re.IGNORECASE)

        if strain_id == '':
            strain_id = 'NA'
            self.nonames.append((doc['accession'], doc['original_strain']))
        doc['strain'] = self.build_canonical_name(sero, country, strain_id, year) # New strain name = DENV1234/COUNTRY/STRAIN_ID/YEAR
        doc['sero'] = sero                                                   # Update serotype in standard format

    def pull_metadata(self, doc, **kwargs): # Called within fix_strain above
        if doc['country'] != None:
            country = doc['country']
        else:
            country = 'NA'
        try:
            division = doc['division']
        except:
            division = 'NA'
        if doc['location'] != None:
            location = doc['location']
        else:
            location = 'NA'
        try:
            year = doc['collection_date'].split('-')[0]
        except AttributeError:
            year = 'NA'
        try:
            month = str(int(doc['collection_date'].split('-')[1]))
        except:
            month = 'NA'
        try:
            day = str(int(doc['collection_date'].split('-')[2]))
        except:
            day = 'NA'
        try:
            sero = 'DENV'+str(int(doc['serotype'][-1])) # Check for serotype annotation like 'dengue virus 1234'
        except TypeError:
            sero = 'DENV'
        try:
            country_code = self.country_to_code[country]
        except KeyError:
            country_code = 'NA'
        return [country_code, country, division, location, sero, year, month, day]

    def build_canonical_name(self, sero, country, strain_id, year):
        return ('%s/%s/%s/%s'%(sero, country, strain_id, year)).strip().upper()

    def separate_viruses_sequences(self, data, **kwargs):
        viruses = []
        sequences = []
        for record in data:
            v = {k: v for k,v in record.items() if k in virus_attribs} # defined in __main__ below
            s = {k: v for k,v in record.items() if k in sequence_attribs}
            v = self.add_virus_fields(v, **kwargs) # add attributes specified at command line, and universal fields like 'number of sequences'
            s = self.add_sequence_fields(s, **kwargs)
            sequences.append(s)
            viruses.append(v)
        return (viruses, sequences)


if __name__=="__main__":
    args = parser.parse_args() # parser is an argparse object initiated in parse.py
    virus_attribs = ['strain', 'original_strain', 'virus', 'serotype','collection_date', 'region', 'country', 'division', 'location'] # define fields in fasta headers that you want used in parse.py > parse > parse_fasta_file ---> (viruses, sequences)
    sequence_attribs = ['accession', 'strain', 'original_strain', 'virus', 'serotype',  'locus', 'sequence', 'authors', 'PMID', 'source', 'gene_list']
    if args.fname == None:
        setattr(args, 'fname', 'results.tbl')
        setattr(args, 'ftype', 'tsv')
    if args.virus == None:
        setattr(args, 'virus', 'dengue')
    if args.database == None:
        setattr(args, 'database', 'vdb')

    setattr(args, 'virus_attribs', virus_attribs)
    setattr(args, 'sequence_attribs', sequence_attribs)

    connVDB = dengue_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
