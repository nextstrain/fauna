import os, argparse, sys
import rethinkdb as r
sys.path.append('')  # need to import from base
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('-edb', '--export_database', default='vdb', help="database to export from")
parser.add_argument('-ev', '--export_virus', default='zika', help="virus table to export")
parser.add_argument('-idb', '--import_database', default='vdb', help="database to import into")
parser.add_argument('-iv', '--import_virus', default='zika', help="virus table to import to")
parser.add_argument('--sync_to_local', default=False, action="store_true", help="sync external database to local database")
parser.add_argument('--sync_from_local', default=False, action="store_true", help="sync local database to external database")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class sync(object):
    def __init__(self, rethink_host=None, auth_key=None, **kwargs):
        if rethink_host is None:
            try:
                self.rethink_host = os.environ['RETHINK_HOST']
            except:
                raise Exception("Missing rethink host")
        elif rethink_host == "localhost":
            raise Exception("Need the rethinkhost that is not local")
        else:
            self.rethink_host = rethink_host
        if auth_key is not None:
            self.auth_key = auth_key
        else:
            try:
                self.auth_key = os.environ['RETHINK_AUTH_KEY']
            except:
                raise Exception("Missing rethink auth_key")
        self.rethink_interact = rethink_interact()

    def sync_from_local(self, export_virus, import_virus, **kwargs):
        '''
        Sync the local rethinkdb instance to an external rethinkdb instance
        Export documents in local database to external database
        '''
        del kwargs['rethink_host']
        del kwargs['auth_key']
        kwargs['export_table'] = export_virus
        kwargs['import_table'] = import_virus
        self.rethink_interact.sync_from_local(self.rethink_host, self.auth_key, pkey='strain', **kwargs)

    def sync_to_local(self, export_virus, import_virus, **kwargs):
        '''
        Sync the local rethinkdb instance to an external rethinkdb instance
        Export documents in external database to local database
        '''
        del kwargs['rethink_host']
        del kwargs['auth_key']
        kwargs['export_table'] = export_virus
        kwargs['import_table'] = import_virus
        self.rethink_interact.sync_to_local(self.rethink_host, self.auth_key, pkey='strain', **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = sync(**args.__dict__)
    if args.sync_to_local:
        connVDB.sync_to_local(**args.__dict__)
    elif args.sync_from_local:
        connVDB.sync_from_local(**args.__dict__)
