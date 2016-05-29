import os, argparse, sys
import rethinkdb as r
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('--from_table', help='database.table to make copy of documents')
parser.add_argument('--to_table', help='database.table to append to')
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")

class append(object):
    def __init__(self, **kwargs):
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_interact = rethink_interact()

    def append(self, **kwargs):
        '''
        Append documents in database.table to another database.table"
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        self.rethink_interact.append(pkey='strain', **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = append(**args.__dict__)
    connVDB.append(**args.__dict__)

