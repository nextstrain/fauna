import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload
from vdb_upload import parser

class Flu_vdb_upload(vdb_upload):

    def __init__(self,  fasta_fields, fasta_fname, database, virus, source, overwrite, locus=None, vsubtype=None, authors=None, path=None, auth_key=None):
        '''
        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        '''
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences', 'virus', 'subtype']
        self.virus_optional_fields = ['division', 'location']
        vdb_upload.__init__(self, fasta_fields, fasta_fname, database, virus, source, overwrite, locus, vsubtype, authors, path, auth_key)

    def format(self):
        '''
        format virus information in preparation to upload to database table
        '''
        print('Formatting for upload')
        self.define_countries()
        self.define_regions()
        for virus in self.viruses:
            self.format_sequence_schema(virus)
            self.format_date(virus)
            self.format_country(virus)
            self.format_region(virus)
            self.format_place(virus)

        # filter out viruses without correct dating format or without region specified
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']), self.viruses)
        self.viruses = filter(lambda v: v['region'] != 'Unknown', self.viruses)
        self.check_all_attributes()

    def define_countries(self):
        '''
        open synonym to country dictionary
        Location is to the level of country of administrative division when available
        '''
        reader = csv.DictReader(open("source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
        self.label_to_country = {}
        for line in reader:
            self.label_to_country[line['label'].lower()] = line['country']

    def format_country(self, v):
        '''
        Label viruses with country based on strain name
        '''
        if "country" not in v:
            v['country'] = 'Unknown'
            try:
                label = re.match(r'^[AB]/([^/]+)/', v['strain']).group(1).lower()						# check first for whole geo match
                if label in self.label_to_country:
                    v['country'] = self.label_to_country[label]
                else:
                    label = re.match(r'^[AB]/([^\-^\/]+)[\-\/]', v['strain']).group(1).lower()			# check for partial geo match
                if label in self.label_to_country:
                    v['country'] = self.label_to_country[label]
                else:
                    label = re.match(r'^[AB]/([A-Z][a-z]+)[A-Z0-9]', v['strain']).group(1).lower()			# check for partial geo match
                if label in self.label_to_country:
                    v['country'] = self.label_to_country[label]
                if v['country'] == 'Unknown':
                    print("couldn't parse location for", v['strain'])
            except:
                print("couldn't parse location for", v['strain'])

if __name__=="__main__":

    args = parser.parse_args()
    fasta_fields = {0:'strain', 1:'accession', 2:'date', 5:'date', 8:'locus'}
    run = Flu_vdb_upload(fasta_fields, fasta_fname=args.fname, database=args.database, virus=args.virus, overwrite=args.overwrite, source=args.source, locus=args.locus, vsubtype=args.subtype, authors=args.authors, path=args.path)
    run.upload()