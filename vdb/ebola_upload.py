import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

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
    args = parser.parse_args()
    virus_fasta_fields = {2:'strain', 3: 'authors', 4:'country', 5:'division', 8:'collection_date'}
    sequence_fasta_fields = {1:'accession', 2:'strain'}
    # 0    1        2                   3            4   5           6 7           8
    #>EBOV|KU296357|2208_C2_18642R_EMLH|EMLab/UoC/WT|SLE|WesternArea|?|Ion_torrent|2015-07-08
    # 0    1                  2                  3        4   5       6 7      8
    #>EBOV|EM_COY_2015_014098|EM_COY_2015_014098|EMLab-RT|GUI|Conakry|?|MinION|2015-03-26
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = ebola_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
