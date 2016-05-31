import os, argparse, sys
import rethinkdb as r
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('--local_table', help='local database.table to sync')
parser.add_argument('--remote_table', help='remote database.table to sync')
parser.add_argument('--push', default=False, action="store_true", help="push local database documents to remote database")
parser.add_argument('--pull', default=False, action="store_true", help="pull remote database documents to local database")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class sync(object):
    def __init__(self, **kwargs):
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_interact = rethink_interact()

    def push(self, **kwargs):
        '''
        Sync the local rethinkdb instance to an external rethinkdb instance
        Export documents in local database to external database
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        self.rethink_interact.push(**kwargs)

    def pull(self, **kwargs):
        '''
        Sync the local rethinkdb instance to an external rethinkdb instance
        Export documents in external database to local database
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        self.rethink_interact.pull(**kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = sync(**args.__dict__)
    if args.pull:
        connVDB.pull(**args.__dict__)
    elif args.push:
        connVDB.push(**args.__dict__)
