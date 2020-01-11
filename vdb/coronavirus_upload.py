import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class coronavirus_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country', 7:'virus_species'}
    sequence_fasta_fields = {0:'accession', 1:'strain'}
    # 0        1                2  3  4       5  6  7
    #>CS362782|UNKNOWN_CS362782|NA|NA|Unknown|NA|NA|Severe_acute_respiratory_syndrome_related_coronavirus
    # 0        1              2  3          4     5            6  7
    #>KF186564|Al_Hasa_4_2013|NA|2013_05_01|Human|Saudi_Arabia|NA|Middle_East_respiratory_syndrome_coronavirus
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = coronavirus_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
