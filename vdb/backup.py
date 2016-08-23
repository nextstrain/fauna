import os, argparse, sys, time, datetime
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database from which backup file is/was made")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('--continuous_backup', default=False, action="store_true",  help="continuously backup database to S3")
parser.add_argument('--backup_s3', default=False, action="store_true", help="backup database to s3")
parser.add_argument('--s3_bucket_name', default='vdb-backups', help="name of s3 bucket where backups stored")
parser.add_argument('--backup_local', default=False, action="store_true", help="backup database locally")
parser.add_argument('--path', default='backups', help="path to store backup of database locally")
parser.add_argument('--days_to_expiration', default=50, help="number of days before deleting backup files")

class backup(object):
    def __init__(self, database, **kwargs):
        self.upload_hour = 3
        self.database = database.lower()
        if 'database' in kwargs:
            kwargs['database'] = self.database

        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_interact = rethink_interact()

    def backup(self, backup_s3, backup_local, **kwargs):
        '''
        Backup every table in specified database
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        if backup_s3:
            self.rethink_interact.backup_s3(**kwargs)
        elif backup_local:
            self.rethink_interact.backup_local(**kwargs)
        else:
            raise Exception("Please specify whether to backup to s3 or locally")

    def continuous_backup(self, **kwargs):
        '''
        Continuously run backup script, uploading tar files to s3 every 24 hours
        '''
        print("Will backup " + self.database + " every 24 hours at hour " + str(self.upload_hour))
        while True:
            if self.time(self.upload_hour):
                self.backup(**kwargs)
                print("Waiting for the next upload hour: " + str(self.upload_hour))
            time.sleep(3600)

    def time(self, upload_hour):
        return int(datetime.datetime.now().hour) == upload_hour

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = backup(**args.__dict__)
    if args.continuous_backup:
        connVDB.continuous_backup(**args.__dict__)
    else:
        connVDB.backup(**args.__dict__)