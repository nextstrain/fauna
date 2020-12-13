import os, argparse, sys
from rethinkdb import r
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact
from vdb.append import parser

class append(object):
    def __init__(self, **kwargs):
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_interact = rethink_interact()

    def append(self, virus, from_database='vdb', to_database='test_vdb', **kwargs):
        '''
        Append documents in database.table to another database.table"
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        kwargs['key'] = 'index'
        virus = virus.lower()
        from_table, to_table = from_database+"."+virus, to_database+"."+virus
        self.rethink_interact.append(from_table=from_table, to_table=to_table, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = append(**args.__dict__)
    connVDB.append(**args.__dict__)
