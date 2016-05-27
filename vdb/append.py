import os, argparse, sys
import rethinkdb as r
sys.path.append('')  # need to import from base
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('--from_table', help='database.table to make copy of documents')
parser.add_argument('--to_table', help='database.table to append to')
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")

class append(object):
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

