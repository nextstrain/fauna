import os, argparse, sys, time, datetime
sys.path.append('')  # need to import from base
print(sys.path)
from vdb.backup import backup

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

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = backup(**args.__dict__)
    if args.continuous_backup:
        connVDB.continuous_backup(**args.__dict__)
    else:
        connVDB.backup(**args.__dict__)