import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class zika_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def fix_name(self, name):
        tmp_name = name
        tmp_name = tmp_name.replace('Human', '').replace('human', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
        tmp_name = tmp_name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
        tmp_name = tmp_name.replace('Zika_virus', '').replace('Zika virus', '').replace('Zika', '').replace('ZIKV', '')
        tmp_name = tmp_name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
        tmp_name = re.sub('^/', '', tmp_name)
        try:
            tmp_name = 'V' + str(int(tmp_name))
        except:
            pass
        return tmp_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    args = parser.parse_args()
    virus_fasta_fields = {2:'strain', 4:'collection_date', 5: 'host', 6:'country'}
    sequence_fasta_fields = {0:'accession', 2:'strain'}
    # 0        1          2      3  4          5     6
    #>KU501216|Zika_virus|103344|NA|2015_12_01|Human|Guatemala
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = zika_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
