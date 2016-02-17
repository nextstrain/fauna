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

        # fields that are needed to upload
        self.upload_fields = ['strain', 'date', 'country', 'gene_region', 'seq']
        self.viruses = self.parse_fasta(self.fasta_file)
        self.strain_name_location = True
        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key="KeHiybPoX8BM6aAhfYqy").repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception
        #create virus table
        #r.db(self.database).table_create("H3N2").run()

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
                v['seq'] = str(record.seq)
                viruses.append(v)
            handle.close()
        return viruses

    def upload(self):
        self.format()
        #self.upload_documents()

    def format(self):
        self.format_date()
        # filter out viruses without correct dating format, must at least have the year specified
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']) != None, self.viruses)

        if self.strain_name_location:
            self.filter_geo()
        '''
        for virus in self.viruses:
            for key, value in virus.items():
                print(key)
                print(value)
            #print(virus['country'] + "\n")
            #print(virus['region'] + "\n")
            print(virus.keys())
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
        Check that each virus has at least the strain name, date, country, sequence region, and sequence.
        Filter out those that don't
        :return:
        '''
        self.viruses = filter(lambda v: all(atr in v for atr in self.upload_fields), self.viruses)

    def upload_documents(self):
        '''
        Insert viruses into collection
        :return:
        '''
        for virus in self.viruses:
            try:
                r.table(self.virus_type).insert(virus).run()
                print("Successfully inserted viruses into database")
            except:
                print("Failed to upload viruses to database")
                raise Exception


fasta_fields = {0:'strain', 5:'date', 8:"gene_region"}
fasta_file = 'data/H3N2_gisaid_epiflu_sequence.fasta'
upload = vdb_upload(fasta_fields, fasta_file)
upload.upload()
