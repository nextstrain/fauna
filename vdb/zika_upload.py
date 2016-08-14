import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class zika_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.strain_fix_fname = "source-data/zika_strain_name_fix.tsv"

    def fix_name(self, doc):
        original_name = doc['strain']
        tmp_name = self.replace_strain_name(original_name, self.fix_whole_name)
        tmp_name = tmp_name.replace('Human', '').replace('human', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
        tmp_name = tmp_name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
        tmp_name = tmp_name.replace('Zika_virus', '').replace('Zikavirus', '').replace('Zika virus', '').replace('Zika', '').replace('ZIKV', '')
        tmp_name = tmp_name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
        tmp_name = re.sub('^/', '', tmp_name)
        try:
            tmp_name = 'V' + str(int(tmp_name))
        except:
            pass
        tmp_name = self.replace_strain_name(tmp_name, self.fix_whole_name)
        return tmp_name, original_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    args = parser.parse_args()
    virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country'}
    sequence_fasta_fields = {0:'accession', 1:'strain'}
    # 0        1      2  3          4     5         6  7
    #>KU501216|103344|NA|2015_12_01|Human|Guatemala|NA|Zika_virus
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = zika_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
