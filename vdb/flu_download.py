import os,datetime
from download import download
from download import get_parser

class flu_download(download):

    def __init__(self, **kwargs):
        download.__init__(self, **kwargs)
        self.virus_specific_fasta_fields = ['vtype', 'subtype', 'lineage']


if __name__=="__main__":
    parser = get_parser()
    parser.add_argument('--subtype', default='any', help='viral lineage')
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.table + '_' + ((args.subtype+'_') if args.subtype!='any' else '') + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connfluVDB = flu_download(**args.__dict__)
    connfluVDB.download(output=False, **args.__dict__)
    if args.subtype!='any':
        connfluVDB.viruses = filter(lambda x:x['subtype']==args.subtype, connfluVDB.viruses)
        print('Documents after subtype filtering:',len(connfluVDB.viruses))
    connfluVDB.output(**args.__dict__)
