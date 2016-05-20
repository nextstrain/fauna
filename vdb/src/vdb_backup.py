import os, argparse, datetime, shutil, subprocess, re, time
import rethinkdb as r
import boto3

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to make backup of")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--backup', default=False, action="store_true", help="backup database to S3")
parser.add_argument('--continuous', default=False, action="store_true",  help="continuously backup database to S3")
parser.add_argument('--restore', default=False, action="store_true", help="restore database to previous version")
parser.add_argument('--restore_table', default=None, help="table to restore")
parser.add_argument('--restore_date', default=None, help="date to restore table to")

class vdb_backup(object):
    def __init__(self, database, rethink_host=None, auth_key=None, **kwargs):
        self.days_to_expiration = 40
        self.upload_hour = 3
        self.database = database.lower()

        if rethink_host is None:
            try:
                self.rethink_host = os.environ['RETHINK_HOST']
            except:
                raise Exception("Missing rethink host")
        else:
            self.rethink_host = rethink_host
        if auth_key is None:
            try:
                self.auth_key = os.environ['RETHINK_AUTH_KEY']
            except:
                raise Exception("Missing rethink auth_key")
        else:
            self.auth_key = auth_key
        self.connect_rethink()

        if 'AWS_ACCESS_KEY_ID' not in os.environ:
            raise Exception("AWS_ACCESS_KEY_ID")
        if 'AWS_SECRET_ACCESS_KEY' not in os.environ:
            raise Exception("AWS_SECRET_ACCESS_KEY")

        self.connect_rethink()
        self.connect_S3()

    def connect_rethink(self):
        '''
        Connect to rethink database,
        '''
        try:
            r.connect(host=self.rethink_host, port=28015, db=self.database, auth_key=self.auth_key).repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

        self.tables = r.db(self.database).table_list().run()

    def connect_S3(self):
        '''
        Connect to AWS S3, assumes AWS_SECRET_ACCESS_KEY and AWS_ACCESS_KEY_ID in os.environ
        :return:
        '''
        s3 = boto3.resource('s3')
        self.bucket = s3.Bucket('vdb-backups')

    def get_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))

    def continuous_backup(self):
        '''
        Continuously run backup script, uploading tar files to s3 every 24 hours
        '''
        while True:
            if self.time():
                self.backup()
                print("Waiting for the next upload hour: " + str(self.upload_hour))
            time.sleep(3600)

    def time(self):
        return int(datetime.datetime.now().hour) == self.upload_hour

    def backup(self):
        '''
        make copy of input database table and import into another database table
        '''
        print("Backing up " + self.database + " on " + self.get_date())
        if not os.path.isdir('temp'):
            os.makedirs('temp')
        for table in self.tables:
            dump_file = self.get_date() + '_' + self.database + '_' + table + '.tar.gz'
            fname = 'temp/' + dump_file
            self.dump(table, fname)
            self.bucket.upload_file(fname, dump_file)
        shutil.rmtree('temp')
        print("Successfully backed up")

        print("Checking for files within " + str(self.days_to_expiration) + " days to delete from s3")
        for object in self.bucket.objects.all():
            fname = str(object.key)
            date = re.match(r'^[^_]*', str(fname)).group(0)
            if self.expired(date):
                self.bucket.delete_objects(Delete={'Objects': [{'Key': fname}]})
                print("Deleting " + fname)

    def dump(self, table, dump_file):
        '''
        backup self.export_database table to self.backup_file
        :return:
        '''
        command = ['rethinkdb', 'dump', '-c', self.rethink_host, '-a', self.auth_key, '-e', self.database + '.' + table, '-f', dump_file]
        try:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
        except:
            raise Exception("Couldn't dump tar file, make sure " + dump_file + " doesn't exist")

    def expired(self, fdate):
        '''
        determine if the file was created past a certain number of days
        '''
        fdate = datetime.datetime.strptime(fdate, '%Y-%m-%d')
        cdate = datetime.datetime.strptime(self.get_date(), '%Y-%m-%d')
        days_since = (cdate-fdate).days
        return days_since >= self.days_to_expiration

    def restore(self, restore_table, restore_date, **kwargs):
        print("Restoring database " + self.database + " and table " + restore_table.lower() + " on to date " + restore_date)
        restore_location = self.database+'.'+restore_table
        restore_file = restore_date + "_" + self.database + "_" + restore_table.lower() + '.tar.gz'
        command = ['rethinkdb', 'restore', restore_file, '-i', restore_location]
        try:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
        except:
            raise Exception("Couldn't restore " + restore_location + " with file " + restore_file)


if __name__=="__main__":
    args = parser.parse_args()
    connVDB = vdb_backup(**args.__dict__)
    if args.continuous:
        connVDB.continuous_backup()
    elif args.restore:
        connVDB.restore(**args.__dict__)
    elif args.backup:
        connVDB.backup()