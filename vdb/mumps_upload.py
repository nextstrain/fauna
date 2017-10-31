import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class mumps_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.strain_fix_fname = "source-data/mumps_strain_name_fix.tsv"
        self.location_fix_fname = "source-data/mumps_location_fix.tsv"
        self.date_fix_fname = "source-data/mumps_date_fix.tsv"

    def fix_name(self, name):
        original_name = name
        name = self.replace_strain_name(original_name, self.fix_whole_name)
        name = name.replace('MuV/', '').replace('MuVi/', '').replace('MuVs/','')
        name = re.sub(r'[_ ]?\[([A-Z])\]$', r'/\1', name)
        name = re.sub(r'\(([A-Z])\)$', r'/\1', name)
        name = re.sub(r'_([A-Z])_$', r'/\1', name)
        name = re.sub(r'[ ;]', r'_', name)
        name = re.sub(r'//', r'/', name)
        name = self.replace_strain_name(name, self.fix_whole_name)
        return name, original_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    args = parser.parse_args()
    virus_fasta_fields = {1:'strain', 2:'collection_date', 3: 'host', 4:'country', 5:'division'}
    sequence_fasta_fields = {0:'accession', 1:'strain'}
    # 0                          1                                2          3     4   5             6
    #>Massachusetts_outbreak_123|MuVs/Massachusetts.USA/50.16/1/G|2016-12-14|Human|USA|Massachusetts|G
    # 0       1                                  2         3     4      5                6
    #>BCCDC90|MuVs/BritishColumbia.CAN/34.16/1/G|2016-8-19|human|canada|british columbia|G
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = mumps_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
