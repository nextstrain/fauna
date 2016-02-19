import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO



class vdb_upload(object):

    def __init__(self, fasta_fields, fasta_file):
        '''

        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        :param fasta_file: File to upload viruses from
        :return: list of viruses
        '''
        print("Uploading Viruses to VDB")

        self.database = "test"
        self.virus_type = "H3N2"
        self.fasta_fields = fasta_fields
        self.fasta_file = fasta_file
        self.virus_source = "GISAID"

        # fields that are needed to upload
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences']
        self.virus_optional_fields = ['division', 'location']
        self.sequence_upload_fields = ['source', 'locus', 'sequence']
        self.sequence_optional_fields = ['accession']  # ex. if from virological.org or not in a database

        self.viruses = self.parse_fasta(self.fasta_file)
        self.strain_name_location = True

        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key="KeHiybPoX8BM6aAhfYqy").repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception
        #create virus table with primary key 'strain'
        try:
            existing_tables = r.db(self.database).table_list().run()
        except:
            print("Couldn't get list of tables")
            raise Exception
        if self.virus_type not in existing_tables:
            try:
                r.db(self.database).table_create(self.virus_type, primary_key='strain').run()
            except:
                print("Couldn't create new table for " + self.virus_type)

    #

    def print_virus_info(self, virus):

        print("----------")
        for key, value in virus.items():
            print(key + " : " + str(value))
        '''
        #print(virus['country'] + "\n")
        #print(virus['region'] + "\n")
        print(virus.keys())
        '''


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
                v['sequence'] = str(record.seq)
                v['virus'] = self.virus_type
                viruses.append(v)
            handle.close()
            print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return viruses

    def format_sequence_schema(self):
        '''
        move sequence information into nested 'sequences' field
        :return:
        '''
        sequence_fields = self.sequence_upload_fields + self.sequence_optional_fields
        for v in self.viruses:
            v['sequences'] = [{}]
            v['sequences'][0]['source'] = self.virus_source
            for field in sequence_fields:
                if field in v.keys():
                    v['sequences'][0][field] = v[field]
                    del v[field]

    def upload(self):
        self.format()
        self.upload_documents()

    def format(self):
        self.format_sequence_schema()

        self.format_date()

        if self.strain_name_location:
            self.filter_geo()
        '''
        for virus in self.viruses:
            self.print_virus_info(virus)
        '''
        self.check_all_attributes()

    def format_date(self):
        '''
        Format viruses date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        '''

        for virus in self.viruses:
            # ex. 2002 (Month and day unknown)
            if re.match(r'\d\d\d\d-\d\d-\d\d', virus['date']):
                pass
            elif re.match(r'\d\d\d\d\s\(Month\sand\sday\sunknown\)', virus['date']):
                virus['date'] = virus['date'][0:4] + "-XX-XX"
            # ex. 2009-06 (Day unknown)
            elif re.match(r'\d\d\d\d-\d\d\s\(Day\sunknown\)', virus['date']):
                virus['date'] = virus['date'][0:7] + "-XX"
            else:
                print("Couldn't reformat this date: " + virus['date'])


        # filter out viruses without correct dating format, must at least have the year specified
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']) != None, self.viruses)

    def filter_geo(self, prune = True):
        '''
        Label viruses with country based on strain name
        Label viruses with region based on country, if prune then filter out viruses without region
        '''
        """Location is to the level of country of administrative division when available"""
        reader = csv.DictReader(open("source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
        label_to_country = {}
        for line in reader:
            label_to_country[line['label'].lower()] = line['country']
        for v in self.viruses:
            if "country" not in v:
                v['country'] = 'Unknown'
                try:
                    label = re.match(r'^[AB]/([^/]+)/', v['strain']).group(1).lower()						# check first for whole geo match
                    if label in label_to_country:
                        v['country'] = label_to_country[label]
                    else:
                        label = re.match(r'^[AB]/([^\-^\/]+)[\-\/]', v['strain']).group(1).lower()			# check for partial geo match
                    if label in label_to_country:
                        v['country'] = label_to_country[label]
                    else:
                        label = re.match(r'^[AB]/([A-Z][a-z]+)[A-Z0-9]', v['strain']).group(1).lower()			# check for partial geo match
                    if label in label_to_country:
                        v['country'] = label_to_country[label]
                    if v['country'] == 'Unknown':
                        print("couldn't parse location for", v['strain'])
                except:
                    print("couldn't parse location for", v['strain'])

        reader = csv.DictReader(open("source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
        country_to_region = {}
        for line in reader:
            country_to_region[line['country']] = line['region']
        for v in self.viruses:
            v['region'] = 'Unknown'
            if v['country'] in country_to_region:
                v['region'] = country_to_region[v['country']]
            if v['country'] != 'Unknown' and v['region'] == 'Unknown':
                print("couldn't parse region for", v['strain'], "country:", v["country"])

        if prune:
            self.viruses = filter(lambda v: v['region'] != 'Unknown', self.viruses)

    def check_all_attributes(self):
        '''
        Assigns 'None' to optional attributes that are missing
        Checks that each virus has upload attributes, filters out viruses that don't.
        :return:
        '''
        self.assign_none()
        self.viruses = filter(lambda v: self.check_upload_attributes(v), self.viruses)

    def assign_none(self):
        for virus in self.viruses:
            for atr in self.virus_optional_fields:
                if atr not in virus:
                    virus[atr] = None
            for atr in self.sequence_optional_fields:
                if atr not in virus['sequences'][0]:
                    virus['sequences'][0][atr] = None

    def check_upload_attributes(self, virus):
        '''

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

            print("Inserting next virus into database: " + virus['strain'])
            try:
                document = r.table(self.virus_type).get(virus['strain']).run()
            except:
                print("Couldn't get document from table")
                print("Strain: " + virus['strain'])
                raise Exception
            if document is None:
                try:
                    r.table(self.virus_type).insert(virus).run()
                    print("Successfully inserted virus into database")
                except:
                    print("Couldn't insert new virus")
                    print("Strain: " + virus['strain'])
            else:
                print(document)
                doc_seqs = document['sequences']
                virus_seq = virus['sequences']
                print(virus_seq[0]['accession'])
                if all(doc_sequence_info['accession'] != virus_seq[0]['accession'] and
                       doc_sequence_info['sequence'] != virus_seq[0]['sequences']
                       for doc_sequence_info in doc_seqs):
                    try:
                        r.table(self.virus_type).get(virus['strain']).update({"sequences": r.row["sequences"].append(
                            virus_seq[0])}).run()
                    except:
                        print("Failed to update virus")
                        print("Strain: " + virus['strain'])
                        raise Exception
                    print("Successfully inserted virus into database")
                else:
                    print("This virus accession and sequence already exists in the database")
                    print("Strain: " + virus['strain'])
                    print("Accession: " + virus_seq[0]['accession'])



fasta_fields = {0:'strain', 1:'accession', 5:'date', 8:"locus"}
fasta_file = 'data/H3N2_gisaid_epiflu_sequence.fasta'
upload = vdb_upload(fasta_fields, fasta_file)
upload.upload()
