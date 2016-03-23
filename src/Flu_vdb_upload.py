import os, re, time, datetime, csv, sys
import numpy as np
import rethinkdb as r
from Bio import SeqIO
from Bio import AlignIO
from vdb_upload import vdb_upload
from vdb_upload import parser

class Flu_vdb_upload(vdb_upload):

    def __init__(self,  fasta_fields, **kwargs):
        '''
        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        '''
        self.virus_upload_fields = ['strain', 'date', 'country', 'sequences', 'virus', 'date_modified', 'subtype']
        self.virus_optional_fields = ['division', 'location']
        vdb_upload.__init__(self, fasta_fields, **kwargs)

        self.patterns = {('A / H3N2', ''): 'H3N2',
                    ('A / H1N1', 'pdm09'): 'H1N1pdm',
                    ('B / H0N0', 'Victoria'): 'Vic',
                    ('B / H0N0', 'Yamagata'): 'Yam',
                    ('A / H1N1', 'seasonal'): 'H1N1',
                    ('A / H7N9', ''): 'H7N9',
                    ('A / H5N1', ''): 'H5N1',
                    ('A / H6N1', ''): 'H6N1',
                    ('A / H5N6', ''): 'H5N6'}
        self.outgroups = {lineage: SeqIO.read('source-data/'+lineage+'_outgroup.gb', 'genbank') for lineage in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']}

    def parse_fasta_file(self, fasta):
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
                v['date_modified'] = self.get_upload_date()
                if 'locus' not in v and self.locus is not None:
                    v['locus'] = self.locus.title()
                if 'authors' not in v and self.authors is not None:
                    v['authors'] = self.authors.title()
                if self.vsubtype is not None:
                    v['subtype'] = self.vsubtype.title()
                else:
                    v['subtype'] = self.determine_lineage(record, v['strain'], (v['subtype'], v['lineage']))
                del v['lineage']
                if 'source' not in v and self.virus_source is not None:
                    v['source'] = self.virus_source.title()
                viruses.append(v)
            handle.close()
            print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return viruses

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
            self.format_place(virus)
            self.format_region(virus)

        # filter out viruses without correct dating format or without region specified
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']), self.viruses)
        self.viruses = filter(lambda v: v['region'] != '?', self.viruses)
        self.viruses = filter(lambda v: v['subtype'] != '?', self.viruses)
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
            v['country'] = '?'
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
                    print("couldn't parse country for", v['strain'])
            except:
                print("couldn't parse country for", v['strain'])

    def determine_lineage(self, seq, strain, tmp_lineage):
        if tmp_lineage in self.patterns:
            return self.patterns[tmp_lineage]
        else:
            scores = []
            for olineage, oseq in self.outgroups.items():
                SeqIO.write([oseq, seq], "temp_in.fasta", "fasta")
                os.system("mafft --auto temp_in.fasta > temp_out.fasta 2>tmp")
                tmp_aln = np.array(AlignIO.read('temp_out.fasta', 'fasta'))
                scores.append((olineage, (tmp_aln[0]==tmp_aln[1]).sum()))
            scores.sort(key = lambda x:x[1], reverse=True)
            if scores[0][1]>0.85*len(seq):
                print(strain, tmp_lineage, len(seq), "\n\t lineage based on similarity:",scores[0][0],"\n\t",scores)
                return scores[0][0]
            else:
                print(strain, tmp_lineage, len(seq), "\n\t other: best scores:", scores[0])
                return '?'

if __name__=="__main__":

    args = parser.parse_args()
    fasta_fields = {0: 'strain', 1: 'accession', 2: 'subtype', 4:'lineage', 5: 'date', 8: 'locus'}
    run = Flu_vdb_upload(fasta_fields, **args.__dict__)
    run.upload()