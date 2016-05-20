import os,datetime
from vdb_download import vdb_download
from vdb_download import parser

class flu_vdb_download(vdb_download):

    def __init__(self, **kwargs):
        vdb_download.__init__(self, **kwargs)
        self.virus_specific_fasta_fields = ['vtype', 'subtype', 'lineage']


if __name__=="__main__":
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connfluVDB = flu_vdb_download(**args.__dict__)
    connfluVDB.download(**args.__dict__)