import os,datetime
from download import download
from download import get_parser
import rethinkdb as r
import time
import re

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
    args.fasta_fields = ['strain', 'accession', 'collection_date', 'region', 'country', 'division', 'location', 'authors']
    if args.virus == None:
        setattr(args, 'virus', 'dengue')
    if args.database == None:
        setattr(args, 'database', 'vdb')
    if args.fstem is None:
        try:
            serotype=args.select[0][-1]
            args.fstem = 'dengue_denv%s'%serotype
        except:
            args.fstem = 'dengue_all'

    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connfluVDB = dengue_download(**args.__dict__)
    connfluVDB.download(**args.__dict__)
