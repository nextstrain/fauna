import os, re, time, datetime, csv, sys
from rethinkdb import r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class zibra_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.strain_fix_fname = "source-data/zika_strain_name_fix.tsv"

    def fix_name(self, name):
        original_name = name
        name = self.replace_strain_name(original_name, self.fix_whole_name)
        name = name.replace('Zika_virus', '').replace('Zikavirus', '').replace('Zika virus', '').replace('Zika', '').replace('ZIKV', '')
        name = name.replace('Human', '').replace('human', '').replace('H.sapiens_wt', '').replace('H.sapiens_tc', '').replace('Hsapiens_tc', '').replace('H.sapiens-tc', '').replace('Homo_sapiens', '').replace('Homo sapiens', '').replace('Hsapiens', '').replace('H.sapiens', '')
        name = name.replace('/Hu/', '')
        name = name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
        name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('__', '_').replace('.', '').replace(',', '')
        name = re.sub('^[\/\_\-]', '', name)
        try:
            name = 'V' + str(int(name))
        except:
            pass
        name = self.replace_strain_name(name, self.fix_whole_name)
        return name, original_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    virus_fasta_fields = {0:'strain', 2:'collection_date', 3:'country', 4: 'division', 5:'location'}
    sequence_fasta_fields = {0:'strain', 1:'accession'}
    # 0                   1                   2          3      4       5
    #>Brazil/ZBRD116/2015|Brazil/ZBRD116/2015|2015-08-28|brazil|alagoas|arapiraca
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = zibra_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
