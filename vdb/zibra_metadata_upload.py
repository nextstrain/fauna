import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class zibra_metadata_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def fix_name(self, name):
        return name

if __name__=="__main__":
    args = parser.parse_args()
#    fasta_fields = {0:'strain', 1:'municipality', 2:'state', 3:'collection_date'}
    fasta_fields = {0:'strain', 1:'location', 2:'division', 3:'date'}
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zibra_metadata_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
