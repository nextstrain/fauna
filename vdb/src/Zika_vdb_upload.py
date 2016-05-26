import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload
from vdb_upload import parser

class Zika_vdb_upload(vdb_upload):

    def __init__(self, **kwargs):
        vdb_upload.__init__(self, **kwargs)
        self.grouping_optional_fields = ['lineage']

    def fix_name(self, name):
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '')
        tmp_name = tmp_name.replace('_Asian', '').replace('_Asia', '').replace('_asian', '').replace('_asia', '')
        try:
            tmp_name = 'V' + str(int(tmp_name))
        except:
            pass
        return tmp_name


if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 2:'strain', 4:'date', 6:'country'}
    # 0        1          2      3  4          5     6
    #>KU501216|Zika_virus|103344|NA|2015_12_01|Human|Guatemala
    setattr(args, 'fasta_fields', fasta_fields)
    if args.path is None:
        args.path = "vdb/data/" + args.virus + "/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = Zika_vdb_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)