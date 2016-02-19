import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload

class H3N2_vdb_upload(vdb_upload):

    def __init__(self, database, virus, source, fasta_fname, fasta_fields):
        vdb_upload.__init__(self, database, virus, source, fasta_fname, fasta_fields)

    def format(self):
        '''
        format virus information in preparation to upload to database table
        :return:
        '''
        self.format_sequence_schema()
        self.format_date()
        self.filter_country()
        self.filter_region()
        self.check_all_attributes()

    def filter_country(self):
        '''
        Label viruses with country based on strain name
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


if __name__=="__main__":

    database = 'test'
    virus = 'H3N2'
    source = 'GISAID'
    fasta_fname = 'H3N2_gisaid_epiflu_sequence.fasta'
    fasta_fields = {0:'strain', 1:'accession', 5:'date', 8:"locus"}


    upload_run = H3N2_vdb_upload(database, virus, source, fasta_fname, fasta_fields)
    upload_run.upload()
