import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload
from vdb_upload import parser

class Flu_vdb_upload(vdb_upload):

    def __init__(self,  fasta_fields, fasta_fname, database, virus, source, locus=None, vsubtype=None, authors=None, path=None, auth_key=None):
        '''
        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        '''
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences', 'virus', 'subtype']
        self.virus_optional_fields = ['division', 'location']
        vdb_upload.__init__(self, fasta_fields, fasta_fname, database, virus, source, locus, vsubtype, authors, path, auth_key)

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
    run = Flu_vdb_upload(fasta_fields, fasta_fname=args.fname, database=args.database, virus=args.virus, source=args.source, locus=args.locus, vsubtype=args.subtype, authors=args.authors, path=args.path)
    run.upload()