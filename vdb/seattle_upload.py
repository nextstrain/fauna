import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class seattle_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    virus_fasta_fields = {0:'strain', 2:'type'}
    sequence_fasta_fields = {0:'strain', 1:'accession', 3:'segment'}
    # 0                                    1                                       2    3
    #>87676fe4-1bf2-4321-ba09-27d81ced2508|87676fe4-1bf2-4321-ba09-27d81ced2508-HA|H3N2|HA
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = seattle_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
