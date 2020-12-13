import os, argparse, sys
from rethinkdb import r
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--virus', help="virus tables")
parser.add_argument('--from_database', help='database.table to make copy of documents')
parser.add_argument('--to_database', help='database.table to append to')
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")

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
        virus = virus.lower()
        viruses_table = virus + "_viruses"
        sequences_table = virus + "_sequences"
        from_table, to_table = from_database+"."+viruses_table, to_database+"."+viruses_table
        self.rethink_interact.append(from_table=from_table, to_table=to_table, **kwargs)
        from_table, to_table = from_database+"."+sequences_table, to_database+"."+sequences_table
        self.rethink_interact.append(from_table=from_table, to_table=to_table, **kwargs)


if __name__=="__main__":
    args = parser.parse_args()
    connVDB = append(**args.__dict__)
    connVDB.append(**args.__dict__)
