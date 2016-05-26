import os,datetime
from download import download
from download import parser

class flu_download(download):

    def __init__(self, **kwargs):
        download.__init__(self, **kwargs)
        self.virus_specific_fasta_fields = ['vtype', 'subtype', 'lineage']


if __name__=="__main__":
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connfluVDB = flu_download(**args.__dict__)
    connfluVDB.download(**args.__dict__)