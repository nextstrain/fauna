import os,datetime
from download import download
from download import get_parser

class seasonal_corona_download(download):
    def __init__(self, **kwargs):
        download.__init__(self, **kwargs)

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    fasta_fields = ['accession', 'strain', 'virus', 'collection_date',
                    'country', 'source', 'sequence_locus', 'subtype', 'type']
    args.fasta_fields = fasta_fields
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connfluVDB = seasonal_corona_download(**args.__dict__)
    connfluVDB.download(**args.__dict__)
