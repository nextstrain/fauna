import os, argparse, sys, time, datetime
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database from which backup file is/was made")
parser.add_argument('-v', '--virus', default=None, help="virus table to restore")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('--s3_bucket_name', default='vdb-backups', help="name of s3 bucket where backups stored")
parser.add_argument('--backup_local', default=False, action="store_true", help="backup database locally")
parser.add_argument('--backup_s3', default=False, action="store_true", help="backup database to s3")
parser.add_argument('--path', default='vdb/backups', help="path to store backup of database locally")
parser.add_argument('--days_to_expiration', default=50, help="number of days before deleting backup files")
parser.add_argument('--restore_date', default=None, help="date to restore table to, format as \'YYYY-MM-DD\'")

class restore(object):
    def __init__(self, database, **kwargs):
        self.upload_hour = 3
        self.database = database.lower()
        if 'database' in kwargs:
            kwargs['database'] = self.database
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_interact = rethink_interact()

    def restore(self, virus, restore_date, **kwargs):
        '''
        restore vdb to previous version
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        virus = virus.lower()
        viruses_table = virus + "_viruses"
        sequences_table = virus + "_sequences"
        self.rethink_interact.restore(restore_table=viruses_table, restore_date=restore_date, **kwargs)
        self.rethink_interact.restore(restore_table=sequences_table, restore_date=restore_date, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = restore(**args.__dict__)
    connVDB.restore(**args.__dict__)