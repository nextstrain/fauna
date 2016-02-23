import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO



class vdb_upload(object):

    def __init__(self, database, virus, source, fasta_fname, fasta_fields, locus=None, path='data/'):
        '''

        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        :param fasta_file: File to upload viruses from
        :return: list of viruses
        '''
        print("Uploading Viruses to VDB")

        self.database = database
        self.virus_type = virus.title()
        self.locus = locus.title()
        self.fasta_fields = fasta_fields
        self.path = path
        self.fasta_file = self.path + self.virus_type + "/" + fasta_fname
        self.virus_source = source

        # fields that are needed to upload
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences', 'virus']
        self.virus_optional_fields = ['division', 'location']
        self.sequence_upload_fields = ['source', 'locus', 'sequence']
        self.sequence_optional_fields = ['accession']  # ex. if from virological.org or not in a database

        self.viruses = self.parse_fasta(self.fasta_file)
        self.strain_name_location = True

        # connect to database
        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key="KeHiybPoX8BM6aAhfYqy").repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

        #create virus table with primary key 'strain'
        existing_tables = r.db(self.database).table_list().run()
        if self.virus_type not in existing_tables:
            r.db(self.database).table_create(self.virus_type, primary_key='strain').run()

        #for formatting region
        try:
            reader = csv.DictReader(open("source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
        except:
            print("Couldn't find geo regions file")
            raise Exception
        self.country_to_region = {}
        for line in reader:
            self.country_to_region[line['country']] = line['region']

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
            print("Tried to find :" + self.fasta_file)
        else:
            for record in SeqIO.parse(handle, "fasta"):
                content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
                v = {key: content[ii] if ii < len(content) else "" for ii, key in self.fasta_fields.items()}
                v['sequence'] = str(record.seq).upper()
                v['virus'] = self.virus_type
                if self.locus is not None:
                    v['locus'] = self.locus
                viruses.append(v)
            handle.close()
            print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return viruses


    def upload(self):
        self.format()
        self.upload_documents()

    def format(self):
        '''
        format virus information in preparation to upload to database table
        :return:
        '''
        for virus in self.viruses:
            self.format_sequence_schema(virus)
            self.format_date(virus)
            self.filter_region(virus)
            self.format_place(virus)
        self.viruses = filter(lambda v: v['region'] != 'Unknown', self.viruses)
        self.check_all_attributes()

    def format_sequence_schema(self, virus):
        '''
        move sequence information into nested 'sequences' field
        :return:
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


        # filter out viruses without correct dating format, must have at least have the year specified
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']), self.viruses)

    def filter_region(self, virus):
        '''
        Label viruses with region based on country, if prune then filter out viruses without region
        :param prune:
        :return:
        '''
        virus['region'] = 'Unknown'
        if virus['country'] in self.country_to_region:
            virus['region'] = self.country_to_region[virus['country']]
        if virus['country'] != 'Unknown' and virus['region'] == 'Unknown':
            print("couldn't parse region for " + virus['strain'] + " country: " + virus["country"])

    def format_place(self, virus):
        '''
        Ensure Camelcase formatting
        :return:
        '''
        location_fields = ['region', 'country', 'division', 'location']
        for field in location_fields:
            virus[field] = virus[field].title()

    def check_all_attributes(self):
        '''
        Assigns 'None' to optional attributes that are missing
        Checks that each virus has upload attributes, filters out viruses that don't.
        :return:
        '''
        self.check_optional_attributes()
        self.viruses = filter(lambda v: self.check_upload_attributes(v), self.viruses)

    def check_optional_attributes(self):
        '''
        Create and assign 'None' to optional attributes that don't exist
        :return:
        '''
        for virus in self.viruses:
            for atr in self.virus_optional_fields:
                if atr not in virus:
                    virus[atr] = None
                elif virus[atr] == '?':
                    virus[atr] = None
            for atr in self.sequence_optional_fields:
                if atr not in virus['sequences'][0]:
                    virus['sequences'][0][atr] = None
                elif virus['sequences'][0][atr] == '?':
                    virus['sequences'][0][atr] = None

    def check_upload_attributes(self, virus):
        '''
        Check that upload attributes are present, print virus information if missing attributes
        :param virus:
        :return: returns true if it has all required upload attributes, else returns false and prints missing attributes
        '''
        missing_attributes = []
        for atr in self.virus_upload_fields:
            if atr not in virus:
                missing_attributes.append(atr)
        for atr in self.sequence_upload_fields:
            if atr not in virus['sequences'][0]:
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
        :return:
        '''
        print("Uploading " + str(len(self.viruses)) + " viruses to the table")
        for virus in self.viruses:
            print("-----------------------")
            print("Inserting next virus into database: " + virus['strain'])
            # Retrieve virus from table to see if it already exists
            document = r.table(self.virus_type).get(virus['strain']).run()

            # Virus doesn't exist in table yet so add it
            if document is None:
                r.table(self.virus_type).insert(virus).run()

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
                r.table(self.virus_type).get(virus['strain']).update({field: virus[field]}).run()

    def update_document_sequence(self, document, virus):
        '''
        Add sequence to sequence list if it's not already there
        '''
        doc_seqs = document['sequences']
        virus_seq = virus['sequences']
        # check that the sequence information isn't already there
        if any(doc_sequence_info['accession'] == virus_seq[0]['accession'] for doc_sequence_info in doc_seqs) or\
                (virus_seq[0]['accession'] == None) and any(doc_sequence_info['sequence'] ==
                virus_seq[0]['sequences'] for doc_sequence_info in doc_seqs):
            print("This virus accession and/or sequence already exists in the database")
            print("Accession: " + str(virus_seq[0]['accession']))
        else:
            r.table(self.virus_type).get(virus['strain']).update({"sequences": r.row["sequences"].append(
                    virus_seq[0])}).run()


if __name__=="__main__":
    fasta_fields = {0:'strain', 1:'accession', 5:'date', 8:"locus"}
    fasta_file = 'data/H3N2_gisaid_epiflu_sequence.fasta'
    upload = vdb_upload("test", "H3N2", "GISAID",fasta_file, fasta_fields)
    upload.upload()
