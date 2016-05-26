import os, argparse, sys
import rethinkdb as r
sys.path.append('')  # need to import from base
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('-edb', '--export_database', default='vdb', help="database to make copy of")
parser.add_argument('-ev', '--export_virus', default='zika', help="virus table to copy")
parser.add_argument('-idb', '--import_database', default='test_vdb', help="database that will be changed")
parser.add_argument('-iv', '--import_virus', default='zika', help="virus table to import to")
parser.add_argument('--move', default=False, action="store_true", help="copy database and table to another database and table")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class move(object):
    def __init__(self, rethink_host=None, auth_key=None, **kwargs):
        if rethink_host is None:
            try:
                self.rethink_host = os.environ['RETHINK_HOST']
            except:
                raise Exception("Missing rethink host")
        else:
            self.rethink_host = rethink_host
        if self.rethink_host == "localhost":
            self.auth_key = None
        elif auth_key is not None:
            self.auth_key = auth_key
        else:
            try:
                self.auth_key = os.environ['RETHINK_AUTH_KEY']
            except:
                raise Exception("Missing rethink auth_key")

        self.rethink_interact = rethink_interact()

    def move(self, export_virus, import_virus, **kwargs):
        '''
        Move documents in database.table to another database.table"
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        kwargs['export_table'] = export_virus
        kwargs['import_table'] = import_virus
        self.rethink_interact.move(pkey='strain', **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = move(**args.__dict__)
    connVDB.move(**args.__dict__)

