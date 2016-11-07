import os,datetime
from download import download
from download import get_parser

class dengue_download(download):
    def __init__(self, **kwargs):
        download.__init__(self, **kwargs)

    def add_selections_command(self, command, selections=[], **kwargs): # Command is an instance of r.table
        '''
        Add selections filter to command
        '''
        if len(selections)>0:
            for sel in selections: # sel like (field, [value1, value2, ...])
                field = sel[0]
                values = sel[1]
                if field == 'gene_list':
                    print "Only downloading documents with one or more of %s in 'gene_list' field."%str(values)
                    command = command.filter(lambda doc: doc[field] in values)
                else:
                    print("Only downloading documents with field \'" + field + "\' equal to one of " + str(values))
                    command = command.filter(lambda doc: r.expr(values).contains(doc[field]))
        return command

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    fasta_fields = ['strain', 'accession']

    args.fasta_fields = fasta_fields
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connfluVDB = dengue_download(**args.__dict__)
    connfluVDB.download(**args.__dict__)
