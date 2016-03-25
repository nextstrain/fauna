import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload
from vdb_upload import parser

class Zika_vdb_upload(vdb_upload):

    def __init__(self, **kwargs):
        vdb_upload.__init__(self, **kwargs)
        self.virus_optional_fields = ['division', 'location']
        self.updateable_virus_fields = ['date', 'country', 'division', 'location', 'virus']

    def canonicalize(self, virus):
        '''
        Canonicalize strain names to consistent format
        '''
        if 'strain' in virus:
            virus['strain'] = self.remove_strings(virus['strain'], ['_Asian', '_Asia'])

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 2:'strain', 4:'date', 6:'country'}
    # 0        1          2      3  4          5     6
    #>KU501216|Zika_virus|103344|NA|2015_12_01|Human|Guatemala
    setattr(args, 'fasta_fields', fasta_fields)
    run = Zika_vdb_upload(**args.__dict__)
    run.upload()