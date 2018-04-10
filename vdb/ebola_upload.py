import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class ebola_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def fix_name(self, name):
        original_name = name
        try:
            name = 'V' + str(int(name))
        except:
            pass
        return name, original_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    #>EBOV|PL7709|KU296426|sierra_leone|PortLoko|2015-06-06
    # 0    1      2        3            4        5
    virus_fasta_fields = {1:'strain', 3:'country', 4:'division', 5:'collection_date'}
    sequence_fasta_fields = {1:'strain', 2:'accession'}
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = ebola_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
