import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
import argparse
from vdb_parse import vdb_parse

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to upload to")
parser.add_argument('-v', '--virus', help="virus table to interact with")
parser.add_argument('--fname', help="input file name")
parser.add_argument('--ftype', default='fasta', help="input file format, default \"fasta\", other is \"genbank\" or \"accession\"")
parser.add_argument('--accessions', default=None, help="comma seperated list of accessions to be uploaded")
parser.add_argument('--source', default=None, help="source of fasta file")
parser.add_argument('--locus', default=None, help="gene or genomic region for sequences")
parser.add_argument('--authors', default=None, help="authors of source of sequences")
parser.add_argument('--subtype', default=None, help="virus subtype, ie \'H3N2\' or \'Yamagata\' for Flu")
parser.add_argument('--private', default=False, action="store_true",  help ="sequences classified as not public")
parser.add_argument('--path', default=None, help="path to fasta file, default is \"data/virus/\"")
parser.add_argument('--host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--overwrite', default=False, action="store_true",  help ="Overwrite fields that are not none")
parser.add_argument('--email', default=None, help="email to access NCBI database via entrez to get virus information")
parser.add_argument('--auto_upload', default=False, action="store_true", help="search genbank for recent sequences")

class vdb_upload(vdb_parse):

    def __init__(self, **kwargs):
        vdb_parse.__init__(self, **kwargs)

        if 'virus' in kwargs:
            self.virus = kwargs['virus'].title()
        if 'database' in kwargs:
            self.database = kwargs['database']
        if 'source' in kwargs:
            self.virus_source = kwargs['source']
        if 'locus' in kwargs:
            self.locus = kwargs['locus']
        if 'subtype' in kwargs:
            self.vsubtype = kwargs['subtype']
        if 'private' in kwargs:
            self.public = not kwargs['private']  #store opposite private
        if 'authors' in kwargs:
            self.authors = kwargs['authors']
        if 'overwrite' in kwargs:
            self.overwrite = kwargs['overwrite']
        if 'fname' in kwargs:
            self.fname = kwargs['fname']
        if 'ftype' in kwargs:
            self.ftype = kwargs['ftype']
        if 'auto_upload' in kwargs:
            self.auto_upload = kwargs['auto_upload']

        if 'path' in kwargs:
            self.path = kwargs['path']
        if self.path is None:
            self.path = "vdb/data/" + self.virus + "/"
        if not os.path.isdir(self.path):
            os.makedirs(self.path)

        # fields that are needed to upload
        self.sequence_upload_fields = ['source', 'locus', 'sequence']
        self.sequence_optional_fields = ['accession']  # ex. if from virological.org or not in a database
        self.citation_optional_fields = ['authors', 'title', 'url']
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences', 'citations', 'virus', 'date_modified', 'public', 'region']
        self.virus_optional_fields = ['division', 'location', 'subtype']
        self.overwritable_virus_fields = ['date', 'country', 'division', 'location', 'virus', 'subtype', 'public']

        if 'host' in kwargs:
            self.host = kwargs['host']
        if 'RETHINK_HOST' in os.environ and self.host is None:
            self.host = os.environ['RETHINK_HOST']
        if self.host is None:
            raise Exception("Missing rethink host")

        if 'auth_key' in kwargs:
            self.auth_key = kwargs['auth_key']
        if 'RETHINK_AUTH_KEY' in os.environ and self.auth_key is None:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")

        self.strains = {}
        self.connect_rethink()

    def connect_rethink(self):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        try:
            r.connect(host=self.host, port=28015, db=self.database, auth_key=self.auth_key).repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

        existing_tables = r.db(self.database).table_list().run()
        if self.virus not in existing_tables:
            raise Exception("No table exists yet for " + self.virus)

        '''
        #create virus table with primary key 'strain'
        existing_tables = r.db(self.database).table_list().run()
        if self.virus not in existing_tables:
            r.db(self.database).table_create(self.virus, primary_key='strain').run()
        '''

    def print_virus_info(self, virus):
        print("----------")
        for key, value in virus.items():
            print(key + " : " + str(value))

    def get_upload_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))

    def upload(self):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to VDB")
        self.parse()
        self.format()
        self.filter()
        self.upload_documents()

    def format(self):
        '''
        format virus information in preparation to upload to database table
        '''
        print('Formatting for upload')
        self.define_regions()
        for virus in self.viruses:
            self.canonicalize(virus)
            self.format_schema(virus)
            self.format_date(virus)
            self.format_place(virus)
            self.format_region(virus)
            self.delete_extra_fields(virus)

    def canonicalize(self, virus):
        '''
        Canonicalize strain names to consistent format
        '''
        if 'strain' in virus:
            virus['strain'] = re.sub(r"\s+", '-', virus['strain'])

    def remove_strings(self, name, strings, replace):
        '''
        Replace each of the strings in name with replace
        '''
        for word in strings:
            name = re.sub(word, replace, name, re.IGNORECASE)
        return name

    def format_schema(self, virus):
        '''
        move sequence information into nested 'sequences' field
        '''
        sequence_fields = self.sequence_upload_fields + self.sequence_optional_fields
        virus['sequences'] = [{}]
        virus['citations'] = [{}]
        for field in sequence_fields:
            if field in virus.keys():
                virus['sequences'][0][field] = virus[field]
                del virus[field]
        for field in self.citation_optional_fields:
            if field in virus.keys():
                virus['citations'][0][field] = virus[field]
                del virus[field]

    def format_date(self, virus):
        '''
        Format viruses date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        Input date could be YYYY_MM_DD, reformat to YYYY-MM-DD
        '''
        # ex. 2002_04_25 to 2002-04-25
        virus['date'] = re.sub(r'_', r'-', virus['date'])
        # ex. 2002 (Month and day unknown)
        if re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', virus['date']):
            pass
        elif re.match(r'\d\d\d\d\s\(Month\sand\sday\sunknown\)', virus['date']):
            virus['date'] = virus['date'][0:4] + "-XX-XX"
        # ex. 2009-06 (Day unknown)
        elif re.match(r'\d\d\d\d-\d\d\s\(Day\sunknown\)', virus['date']):
            virus['date'] = virus['date'][0:7] + "-XX"
        elif re.match(r'\d\d\d\d-\d\d', virus['date']):
            virus['date'] = virus['date'][0:7] + "-XX"
        elif re.match(r'\d\d\d\d', virus['date']):
            virus['date'] = virus['date'][0:4] + "-XX-XX"
        else:
            print("Couldn't reformat this date: " + virus['date'])

    def define_regions(self):
        '''
        open country to region dictionary
        '''
        try:
            reader = csv.DictReader(open("vdb/source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
        except:
            raise Exception("Couldn't find geo regions file")
        self.country_to_region = {}
        for line in reader:
            self.country_to_region[line['country']] = line['region']

    def format_region(self, virus):
        '''
        Label viruses with region based on country, if prune then filter out viruses without region
        '''
        virus['region'] = None
        if virus['country'] in self.country_to_region:
            virus['region'] = self.country_to_region[virus['country']]
        if virus['country'] != '?' and virus['region'] == '?':
            print("couldn't parse region for " + virus['sequences'][0]['accession'] + ", country: " + virus["country"])

    def format_place(self, virus):
        '''
        Ensure Camelcase formatting for geographic information
        '''
        location_fields = ['country', 'division', 'location']
        for field in location_fields:
            if field in virus:
                if len(virus[field].split()) > 1:
                    virus[field] = virus[field].title()
        virus['country'] = virus['country'].replace("_", "").replace(" ", "")

    def filter(self):
        '''
        Filter out viruses without correct dating format or without region specified
        Check  optional and upload attributes
        '''
        print(str(len(self.viruses)) + " viruses before filtering")
        self.check_optional_attributes()
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']), self.viruses)
        self.viruses = filter(lambda v: isinstance(v['public'], (bool)), self.viruses)
        self.viruses = filter(lambda v: v['region'] is not None, self.viruses)
        self.viruses = filter(lambda v: self.check_upload_attributes(v), self.viruses)
        print(str(len(self.viruses)) + " viruses after filtering")

    def check_optional_attributes(self):
        '''
        Reassign unknowns from '?' to 'None'
        Create and assign 'None' to optional attributes that don't exist
        '''
        for virus in self.viruses:
            for key in virus.keys():
                if virus[key] == '?':
                    virus[key] = None
                if isinstance(virus[key], (str)):
                    virus[key] = virus[key].strip
            for key in virus['sequences'][0].keys():
                if virus['sequences'][0][key] == '?':
                    virus['sequences'][0][key] = None
                if isinstance(virus['sequences'][0][key], (str)):
                    virus['sequences'][0][key] = virus['sequences'][0][key].strip
            for key in virus['citations'][0].keys():
                if virus['citations'][0][key] == '?':
                    virus['citations'][0][key] = None
                if isinstance(virus['citations'][0], (str)):
                    virus['citations'][0] = virus['citations'][0].strip

            for atr in self.virus_optional_fields:
                if atr not in virus:
                    virus[atr] = None
            for atr in self.sequence_optional_fields:
                if atr not in virus['sequences'][0]:
                    virus['sequences'][0][atr] = None
            for atr in self.citation_optional_fields:
                if atr not in virus['citations'][0]:
                    virus['citations'][0][atr] = None

    def check_upload_attributes(self, virus):
        '''
        Checks that required upload attributes are present and not equal to None for given virus
        :return: returns true if it has all required upload attributes, else returns false and prints missing attributes
        '''
        missing_attributes = []
        for atr in self.virus_upload_fields:
            if atr in virus and virus[atr] is not None:
                pass
            else:
                missing_attributes.append(atr)
        for atr in self.sequence_upload_fields:
            if atr in virus['sequences'][0] and virus['sequences'][0][atr] is not None:
                pass
            else:
                missing_attributes.append(atr)
        if len(missing_attributes) > 0:
            print("This strain is missing a required attribute and will be removed from upload sequences")
            print(virus['strain'])
            print("Missing attributes: " + str(missing_attributes))
            return False
        else:
            return True

    def delete_extra_fields(self, virus):
        for key in virus.keys():
            if key not in self.virus_upload_fields + self.virus_optional_fields:
                print("Deleting virus info " + key + ": " + virus[key] + " from " + virus['strain'])
                del virus[key]
        for key in virus['sequences'][0].keys():
            if key not in self.sequence_upload_fields + self.sequence_optional_fields:
                print("Deleting sequence info " + key + ": " + virus['sequences'][0][key] + " from " + virus['strain'])
                del virus['sequences'][0][key]
        for key in virus['citations'][0].keys():
            if key not in self.citation_optional_fields:
                print("Deleting citation info " + key + ": " + virus['citations'][0][key] + " from " + virus['strain'])
                del virus['citations'][0][key]

    def upload_documents(self):
        '''
        Insert viruses into collection
        '''
        print("Uploading " + str(len(self.viruses)) + " viruses to the table")
        self.relaxed_strains()
        for virus in self.viruses:            
            # Retrieve virus from table to see if it already exists, try relaxed comparison first
            relaxed_name = self.relax_name(virus['strain'])
            if relaxed_name in self.strains:
                self.strain_name = self.strains[relaxed_name]
            else:
                self.strain_name = virus['strain']
            document = r.table(self.virus).get(self.strain_name).run()
            # Virus doesn't exist in table yet so add it
            if document is None:
                print("Inserting " + virus['strain'] + " into database")
                r.table(self.virus).insert(virus).run()
            # Virus exists in table so just add sequence information and update meta data if needed
            else:
                self.updated = False
                self.update_document_sequence(document, virus)
                self.update_document_meta(document, virus)

    def relaxed_strains(self):
        '''
        Create dictionary from relaxed vdb strain names to actual vdb strain names.
        '''
        strains = {}
        cursor = list(r.db(self.database).table(self.virus).run())
        for doc in cursor:
            strains[self.relax_name(doc['strain'])] = doc['strain']
        self.strains = strains

    def relax_name(self, name):
        '''
        Return the relaxed strain name to compare with
        '''
        name = re.sub(r"-", '', name)
        name = re.sub(r"_", '', name)
        name = re.sub(r"/", '', name)
        return name

    def update_document_meta(self, document, virus):
        '''
        if overwrite is false update doc_virus information only if virus info is different and not null
        if overwrite is true update doc_virus information if virus info is different
        '''
        for field in self.overwritable_virus_fields:
            # update if overwrite and anything
            # update if !overwrite only if document[field] is not none
            if field not in document:
                if field in virus:
                    print("Creating virus field ", field, " assigned to ", virus[field])
                    r.table(self.virus).get(self.strain_name).update({field: virus[field]}).run()
                    document[field] = virus[field]
            elif (self.overwrite and document[field] != virus[field]) or (not self.overwrite and document[field] is None and document[field] != virus[field]):
                if field in virus:
                    print("Updating virus field " + str(field) + ", from \"" + str(document[field]) + "\" to \"" + virus[field]) + "\""
                    r.table(self.virus).get(self.strain_name).update({field: virus[field]}).run()
                    document[field] = virus[field]
                    self.updated = True
        if self.updated:
            r.table(self.virus).get(self.strain_name).update({'date_modified': virus['date_modified']}).run()
            document['date_modified'] = virus['date_modified']

    def update_document_sequence(self, document, virus):
        '''
        Update sequence fields if matching accession or sequence
        Append sequence to sequence list if no matching accession or sequence
        '''
        doc_seqs = document['sequences']
        virus_seq = virus['sequences'][0]
        virus_citation = virus['citations'][0]
        if virus_seq['accession'] != None:
            if all(virus_seq['accession'] != seq_info['accession'] for seq_info in doc_seqs):
                self.append_new_sequence(document, virus_seq, virus_citation)
            else:
                self.update_sequence_citation_field(virus, document, 'accession', self.sequence_upload_fields+self.sequence_optional_fields, self.citation_optional_fields)
        else:
            if all(virus_seq['sequence'] != seq_info['sequence'] for seq_info in doc_seqs):
                self.append_new_sequence(document, virus_seq, virus_citation)
            else:
                self.update_sequence_citation_field(virus, document, 'sequence', self.sequence_upload_fields+self.sequence_optional_fields, self.citation_optional_fields)

        if len(document['sequences']) != len(document['citations']):
            print("Warning: length of list of sequences and citations does not match for " + virus['strain'])

    def append_new_sequence(self, document, virus_seq, virus_citation):
        try:
            r.table(self.virus).get(self.strain_name).update({"sequences": r.row["sequences"].append(virus_seq)}).run()
            r.table(self.virus).get(self.strain_name).update({"citations": r.row["citations"].append(virus_citation)}).run()
        except:
            print("Couldn't append new sequence and citation info for " + document['strain'])
        document['sequences'].append(virus_seq)
        document['citations'].append(virus_citation)
        self.updated = True

    def update_sequence_citation_field(self, virus, document, check_field, sequence_fields, citation_fields):
        '''
        Checks for matching viruses by comparing sequence and accession, updates other attributes if needed
        '''
        updated_sequence = False
        updated_citation = False
        strain = virus['strain']
        doc_seqs = document['sequences']
        virus_seq = virus['sequences'][0]
        index = -1
        for doc_sequence_info in doc_seqs:
            index += 1
            if doc_sequence_info[check_field] == virus_seq[check_field]:
                updated_sequence = self.update_nested_field(strain, sequence_fields, doc_sequence_info, virus_seq)
                virus_citation = virus['citations'][0]
                doc_citation_info = document['citations'][index]
                updated_citation = self.update_nested_field(strain, citation_fields, doc_citation_info, virus_citation)

        if updated_sequence:
            document['sequences'] = doc_seqs
            r.table(self.virus).get(self.strain_name).update({"sequences": doc_seqs}).run()
            self.updated = True
        if updated_citation:
            document['citations'][index] = doc_citation_info
            r.table(self.virus).get(self.strain_name).update({"citations": document['citations']}).run()
            self.updated = True

    def update_nested_field(self, strain, fields, doc_info, virus_info):
        '''
        check for updates to nested sequences and citations fields.
        '''
        updated_sequence = False
        for field in fields:
            if field not in doc_info:
                if field in virus_info:
                    print("Creating field ", field, " assigned to \"", virus_info[field] + "\" for strain: " + strain)
                    doc_info[field] = virus_info[field]
                    updated_sequence = True
            elif (self.overwrite and doc_info[field] != virus_info[field]) or (not self.overwrite and doc_info[field] is None and doc_info[field] != virus_info[field]):
                if field in virus_info:
                    print("Updating field " + str(field) + ", from \"" + str(doc_info[field]) + "\" to \"" + str(virus_info[field]) + "\" for strain: " + strain)
                    doc_info[field] = virus_info[field]
                    updated_sequence = True
        return updated_sequence

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 1:'strain', 2:'date', 4:'country', 5:'division', 6:'location'}
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = vdb_upload(**args.__dict__)
    connVDB.upload()
