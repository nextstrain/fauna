import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-db', '-database', help="database to download from")
parser.add_argument('-v', '-virus', help="virus table to interact with")
parser.add_argument('-fname', help="input file name, default \"fasta\", other is \"json\"")
parser.add_argument('-source', help="source of fasta file")
parser.add_argument('--locus', default=None, help="gene or genomic region for sequences")
parser.add_argument('--subtype', default=None, help="virus subtype, ie \'H3N2\' or \'Yamagata\' for Flu")
parser.add_argument('--path', default=None, help="path to fasta import files, default is \"data/virus/\"")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")


class vdb_upload(object):

    def __init__(self,  fasta_fields, fasta_fname, database, virus, source, locus=None, vsubtype=None, path=None, auth_key=None):
        '''
        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        '''
        print("Uploading Viruses to VDB")

        self.database = database
        self.virus = virus.title()
        if source is not None:
            self.virus_source = source.title()
        if locus is not None:
            self.locus = locus.title()
        self.fasta_fields = fasta_fields

        self.path = path
        if self.path is None:
            self.path = "data/" + self.virus + "/"
        if not os.path.isdir(self.path):
            os.makedirs(self.path)
        self.fasta_fname = fasta_fname
        self.viruses = self.parse_fasta(self.path + self.fasta_fname)

        # fields that are needed to upload
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences', 'virus']
        self.virus_optional_fields = ['division', 'location', ]
        self.sequence_upload_fields = ['source', 'locus', 'sequence']
        self.sequence_optional_fields = ['accession', 'authors']  # ex. if from virological.org or not in a database

        self.auth_key = auth_key
        if 'RETHINK_AUTH_KEY' in os.environ and self.auth_key is None:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")

        self.connect_rethink()

    def connect_rethink(self):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key=self.auth_key).repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

        #create virus table with primary key 'strain'
        existing_tables = r.db(self.database).table_list().run()
        if self.virus not in existing_tables:
            r.db(self.database).table_create(self.virus, primary_key='strain').run()

    def print_virus_info(self, virus):

        print("----------")
        for key, value in virus.items():
            print(key + " : " + str(value))

    def parse_fasta(self, fasta):
        '''
        Parse FASTA file with default header formatting
        :return: list of documents(dictionaries of attributes) to upload
        '''
        viruses = []
        try:
            handle = open(fasta, 'r')
        except IOError:
            print(fasta, "not found")
        else:
            for record in SeqIO.parse(handle, "fasta"):
                content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
                v = {key: content[ii] if ii < len(content) else "" for ii, key in self.fasta_fields.items()}
                v['sequence'] = str(record.seq).upper()
                v['virus'] = self.virus
                if 'locus' not in v and self.locus is not None:
                    v['locus'] = self.locus
                viruses.append(v)
            handle.close()
            print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return viruses

    def upload(self):
        '''
        format virus information, then upload to database
        :return:
        '''
        self.format()
        self.upload_documents()

    def format(self):
        '''
        format virus information in preparation to upload to database table
        '''
        print('Formatting for upload')
        self.define_regions()
        for virus in self.viruses:
            self.format_sequence_schema(virus)
            self.format_date(virus)
            self.format_region(virus)
            self.format_place(virus)

        # filter out viruses without correct dating format or without region specified
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']), self.viruses)
        self.viruses = filter(lambda v: v['region'] != 'Unknown', self.viruses)
        self.check_all_attributes()

    def format_sequence_schema(self, virus):
        '''
        move sequence information into nested 'sequences' field
        '''
        sequence_fields = self.sequence_upload_fields + self.sequence_optional_fields

        virus['sequences'] = [{}]
        virus['sequences'][0]['source'] = self.virus_source.title()
        for field in sequence_fields:
            if field in virus.keys():
                virus['sequences'][0][field] = virus[field]
                del virus[field]

    def format_date(self, virus):
        '''
        Format viruses date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        '''

        # ex. 2002 (Month and day unknown)
        if re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', virus['date']):
            pass
        elif re.match(r'\d\d\d\d\s\(Month\sand\sday\sunknown\)', virus['date']):
            virus['date'] = virus['date'][0:4] + "-XX-XX"
        # ex. 2009-06 (Day unknown)
        elif re.match(r'\d\d\d\d-\d\d\s\(Day\sunknown\)', virus['date']):
            virus['date'] = virus['date'][0:7] + "-XX"
        else:
            print("Couldn't reformat this date: " + virus['date'])

    def define_regions(self):
        '''
        open country to region dictionary
        '''
        try:
            reader = csv.DictReader(open("source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
        except:
            raise Exception("Couldn't find geo regions file")
        self.country_to_region = {}
        for line in reader:
            self.country_to_region[line['country']] = line['region']

    def format_region(self, virus):
        '''
        Label viruses with region based on country, if prune then filter out viruses without region
        '''
        virus['region'] = 'Unknown'
        if virus['country'] in self.country_to_region:
            virus['region'] = self.country_to_region[virus['country']]
        if virus['country'] != 'Unknown' and virus['region'] == 'Unknown':
            print("couldn't parse region for " + virus['strain'] + " country: " + virus["country"])

    def format_place(self, virus):
        '''
        Ensure Camelcase formatting for geographic information
        '''
        location_fields = ['region', 'country', 'division', 'location']
        for field in location_fields:
            if field in virus:
                virus[field] = virus[field].title()

    def check_all_attributes(self):
        '''
        Check that each virus has upload attributes, filters out viruses that don't.
        '''
        self.check_optional_attributes()
        self.viruses = filter(lambda v: self.check_upload_attributes(v), self.viruses)

    def check_optional_attributes(self):
        '''
        Reassign unknowns from '?' to 'None'
        Create and assign 'None' to optional attributes that don't exist
        '''
        for virus in self.viruses:
            for key in virus.keys():
                if virus[key] == '?':
                    virus[key] = None
            for key in virus['sequences'][0].keys():
                if virus['sequences'][0][key] == '?':
                    virus['sequences'][0][key] = None
            for atr in self.virus_optional_fields:
                if atr not in virus:
                    virus[atr] = None
            for atr in self.sequence_optional_fields:
                if atr not in virus['sequences'][0]:
                    virus['sequences'][0][atr] = None

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

    def upload_documents(self):
        '''
        Insert viruses into collection
        '''
        print("Uploading " + str(len(self.viruses)) + " viruses to the table")
        for virus in self.viruses:
            print("-----------------------")
            print("Inserting next virus into database: " + virus['strain'])
            # Retrieve virus from table to see if it already exists
            document = r.table(self.virus).get(virus['strain']).run()
            # Virus doesn't exist in table yet so add it
            if document is None:
                r.table(self.virus).insert(virus).run()
            # Virus exists in table so just add sequence information and update meta data if needed
            else:
                self.update_document_sequence(document, virus)
                self.update_document_meta(document, virus)

    def update_document_meta(self, document, virus):
        '''
        update doc_virus information if virus info is different and not null
        '''
        updateable_virus_fields = ['date', 'country', 'division', 'location', 'virus']
        for field in updateable_virus_fields:
            if virus[field] != None and document[field] != virus[field]:
                print("Updating virus field ", field, ", from ", document[field], " to ", virus[field])
                r.table(self.virus).get(virus['strain']).update({field: virus[field]}).run()

    def update_document_sequence(self, document, virus):
        '''
        Update sequence fields if matching accession or sequence
        Append sequence to sequence list if no matching accession or sequence
        '''
        doc_seqs = document['sequences']
        virus_seq = virus['sequences'][0]
        if virus_seq['accession'] != None:
            self.update_sequence_field(virus, document, 'accession')
        else:
            self.update_sequence_field(virus, document, 'sequence')
        if (virus_seq['accession'] != None and all(virus_seq['accession'] != seq_info['accession'] for seq_info in doc_seqs)) or (virus_seq['accession'] == None and all(virus_seq['sequence'] != seq_info['sequence'] for seq_info in doc_seqs)):
            r.table(self.virus).get(virus['strain']).update({"sequences": r.row["sequences"].append(virus_seq)}).run()
        else:
            print("This virus already exists in the table")

    def update_sequence_field(self, virus, document, check_field):
        doc_seqs = document['sequences']
        virus_seq = virus['sequences'][0]
        for doc_sequence_info in doc_seqs:
            if doc_sequence_info[check_field] == virus_seq[check_field]:
                for field in (self.sequence_upload_fields+self.sequence_optional_fields):
                    if field not in doc_sequence_info:
                        doc_sequence_info[field] = virus_seq[field]
                    if virus_seq[field] != None and doc_sequence_info[field] != virus_seq[field]:
                        print("Updating virus field \"" + field + "\", from \"" + doc_sequence_info[field] + "\" to \"" + virus_seq[field] + "\"")
                        doc_sequence_info[field] = virus_seq[field]
        r.table(self.virus).get(virus['strain']).update({"sequences": doc_seqs}).run()

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 1:'strain', 2:'date', 4:'country', 5:'division', 6:'location'}
    run = vdb_upload(fasta_fields, fasta_fname=args.fname, database=args.db, virus=args.v, source=args.source,
                     locus=args.locus, vsubtype=args.subtype, path=args.path)
    run.upload()

#python vdb_upload.py -database test -v zika -fname zika_virological_02-22-2016.fasta -source Virological --locus Genome