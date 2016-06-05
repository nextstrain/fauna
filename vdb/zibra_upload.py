import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class zibra_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.grouping_optional_fields = ['lineage']

    def fix_name(self, name):
        return name

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'strain'}
    # 0
    #>LACEN_Sample_ID
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zibra_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
