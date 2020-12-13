import os, re, time, datetime, csv, sys
from rethinkdb import r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class dengue_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.strain_fix_fname = "source-data/dengue_strain_name_fix.tsv"

    def fix_name(self, name):
        original_name = name
        name = self.replace_strain_name(original_name, self.fix_whole_name)
        name = name.replace('Dengue_virus', '')
        name = name.replace('Human', '').replace('human', '').replace('H.sapiens_wt', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
        name = name.replace('_URI', '').replace('_SER', '').replace('_PLA', '').replace('_MOS', '').replace('_SAL', '')
        name = name.replace('Aaegypti_wt', 'Aedes_aegypti').replace('Aedessp', 'Aedes_sp')
        name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
        name = re.sub('^[\/\_\-]', '', name)
        try:
            name = 'V' + str(int(name))
        except:
            pass
        name = self.replace_strain_name(name, self.fix_whole_name)
        return name, original_name

    def fix_casing(self, document):
        if 'host' in document and document['host'] is not None:
            document['host'] = self.camelcase_to_snakecase(document['host'])
        if 'serotype' in document and document['serotype'] is not None:
            document['serotype'] = self.camelcase_to_snakecase(document['serotype'])


if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country', 7: 'serotype'}
    sequence_fasta_fields = {0:'accession', 1:'strain'}
    # 0        1     2  3       4     5         6   7
    #>MF033243|03720|NA|2016_01|Human|Singapore|1_V|Dengue_virus_1
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = dengue_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
