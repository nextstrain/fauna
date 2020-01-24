import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class ncov_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.location_fix_fname = "source-data/ncov_location_fix.tsv"

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
    virus_fasta_fields = {0:'strain', 2:'collection_date', 3: 'host', 4:'country', 5:'division', 6:'city', 7:'originating_lab'}
    sequence_fasta_fields = {0:'strain', 1:'accession', 8:'submitting_lab', 9:'authors'}
    # 0                             1              2          3     4     5     6     7                  8                      9
    #>BetaCoV/Wuhan/IVDC-HB-01/2019|EPI_ISL_402119|2019-12-30|Human|China|Hubei|Wuhan|National Institute|National Institute for|Tan et al
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = ncov_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
