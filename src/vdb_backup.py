import os, argparse, datetime, shutil, subprocess, re, time
import rethinkdb as r
import boto3

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='test', help="database to make backup of")
parser.add_argument('--host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--continuous', default=False, action="store_true",  help="Continuously backup database to S3")

class vdb_backup(object):
    def __init__(self, **kwargs):
        self.days_to_expiration = 40
        self.upload_hour = 3

        if 'virus' in kwargs and kwargs['virus'] is not None:
            self.virus = kwargs['virus'].title()
        if 'database' in kwargs:
            self.database = kwargs['database']

        if 'host' in kwargs:
            self.host = kwargs['host']
        if 'RETHINK_HOST' in os.environ and self.host is None:
            self.host = os.environ['RETHINK_HOST']
        if self.host is None:
            raise Exception("Missing rethink host")

        if 'auth_key' in kwargs:
            self.auth_key = kwargs['auth_key']
        if 'RETHINK_AUTH_KEY' in os.environ and self.auth_key is None:
            self.auth_key = os.environ['RETHINK_AUTH_KEY']
        if self.auth_key is None:
            raise Exception("Missing auth_key")

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
            r.connect(host=self.host, port=28015, db=self.database, auth_key=self.auth_key).repl()
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

    def dump(self, table, dump_file):
        '''
        backup self.export_database table to self.backup_file
        :return:
        '''
        command = ['rethinkdb', 'dump', '-c', self.host, '-a', self.auth_key, '-e', self.database + '.' + table, '-f', dump_file]
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

if __name__=="__main__":
    args = parser.parse_args()
    run = vdb_backup(**args.__dict__)
    if args.continuous:
        run.continuous_backup()
    else:
        run.backup()