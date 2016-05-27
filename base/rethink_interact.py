import os, shutil, sys, datetime, re, subprocess
import rethinkdb as r
import boto3
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

class rethink_interact(object):
    def __init__(self, **kwargs):
        self.rethink_io = rethink_io()

    def connect_S3(self, s3_bucket_name, **kwargs):
        '''
        Connect to AWS S3, checks for AWS_SECRET_ACCESS_KEY and AWS_ACCESS_KEY_ID in os.environ
        :return:
        '''
        if 'AWS_ACCESS_KEY_ID' not in os.environ:
            raise Exception("AWS_ACCESS_KEY_ID")
        if 'AWS_SECRET_ACCESS_KEY' not in os.environ:
            raise Exception("AWS_SECRET_ACCESS_KEY")
        if s3_bucket_name is None:
            raise Exception("Please give s3 bucket name to connect to")

        try:
            s3 = boto3.resource('s3')
            return s3.Bucket(s3_bucket_name)
        except:
            raise Exception("Couldn't connect to s3 bucket " + s3_bucket_name)

    def backup_s3(self, database, **kwargs):
        '''
        make backup of every table in database, upload to s3 bucket
        '''
        print("Backing up " + database + " on " + self.rethink_io.get_upload_date())
        if not os.path.isdir('temp'):
            os.makedirs('temp')
        tables = r.db(database).table_list().run()
        bucket = self.connect_S3(**kwargs)
        for table in tables:
            dump_file = self.rethink_io.get_upload_date() + '_' + database + '_' + table + '.tar.gz'
            fname = 'temp/' + dump_file
            self.dump(database, table, fname, **kwargs)
            bucket.upload_file(fname, dump_file)
        shutil.rmtree('temp')
        print("Successfully backed up")
        self.delete_expired_s3_backups(bucket, **kwargs)

    def backup_local(self, database, path='', **kwargs):
        '''
        make backup of every table in database, store locally
        '''
        print("Backing up " + database + " on " + self.rethink_io.get_upload_date() + " to location: " + path)
        if not os.path.isdir(path):
            os.makedirs(path)
        tables = r.db(database).table_list().run()
        for table in tables:
            dump_file = self.rethink_io.get_upload_date() + '_' + database + '_' + table + '.tar.gz'
            self.dump(database=database, dump_table=table, dump_file=dump_file, **kwargs)
            shutil.move(dump_file, path+'/'+dump_file)
            print("Successfully backed up " + table)
        self.delete_expired__local_backups(path=path, **kwargs)

    def dump(self, database, dump_table, dump_file, rethink_host='localhost', auth_key=None, **kwargs):
        '''
        backup self.export_database table to self.backup_file
        :return:
        '''
        if rethink_host!='localhost' and auth_key is not None:
            command = ['rethinkdb', 'dump', '-c', rethink_host, '-a', auth_key, '-e', database + '.' + dump_table, '-f', dump_file]
        else:
            command = ['rethinkdb', 'dump', '-e', database + '.' + dump_table, '-f', dump_file]
        try:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
        except:
            raise Exception("Couldn't dump tar file, make sure " + dump_file + " doesn't exist")

    def restore(self, database, restore_table, restore_date, rethink_host='localhost', auth_key=None, **kwargs):
        '''
        '''
        if restore_table is not None and restore_date is not None:
            print("Restoring database " + database + " and table " + restore_table + " to date " + restore_date)
            restore_location = database+'.'+restore_table
            restore_file = restore_date + "_" + database + "_" + restore_table + '.tar.gz'
            self.get_file(fname=restore_file, **kwargs)
            if rethink_host != 'localhost' and auth_key is not None:
                command = ['rethinkdb', 'restore', restore_file, '-i', restore_location, '-c', rethink_host, '-a', auth_key, '--force']
            else:
                command = ['rethinkdb', 'restore', restore_file, '-i', restore_location, '--force']
            try:
                with open(os.devnull, 'wb') as devnull:
                    subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
                os.remove(restore_file)
            except:
                raise Exception("Couldn't restore " + str(restore_location) + " with file " + str(restore_file))
        else:
            raise Exception("define arguments \'--restore_table\' and \'--restore_date\'")

    def get_file(self, s3, local, fname, path='', **kwargs):
        '''
        Download file from s3 to current directory
        '''
        if s3:
            bucket = self.connect_S3(**kwargs)
            existing_files = [str(object.key) for object in bucket.objects.all()]
            if fname in existing_files:
                bucket.download_file(fname, fname)
            else:
                raise Exception("Specified file could not be found in s3 bucket", fname)
        elif local:
            shutil.copy(path+'/'+fname, fname, )
            if not os.path.isfile(path+'/'+fname):
                raise Exception("Couldn't find backup file at specified location", fname)
        else:
            raise Exception("Please specify whether to restore from a backup stored in s3 or locally")

    def delete_expired_s3_backups(self, bucket, days_to_expiration, **kwargs):
        '''
        delete backup files from s3 that are more than some amount of days old
        '''
        print("Checking for files within " + str(days_to_expiration) + " days to delete from s3")
        for object in bucket.objects.all():
            fname = str(object.key)
            date = re.match(r'^[^_]*', str(fname)).group(0)
            if self.expired(fdate=date, days=days_to_expiration):
                bucket.delete_objects(Delete={'Objects': [{'Key': fname}]})
                print("Deleting " + fname)

    def delete_expired__local_backups(self, path, days_to_expiration, **kwargs):
        '''
        delete backup files stored locally that are more than some amount of days old
        '''
        for fname in os.listdir(path):
            date = re.match(r'^[^_]*', str(fname)).group(0)
            if self.expired(date, days_to_expiration):
                os.remove(path+fname)
                print("Deleting " + fname)

    def expired(self, fdate, days):
        '''
        determine if the file was created past a certain number of days
        '''
        fdate = datetime.datetime.strptime(fdate, '%Y-%m-%d')
        cdate = datetime.datetime.strptime(self.rethink_io.get_upload_date(), '%Y-%m-%d')
        days_since = (cdate-fdate).days
        return days_since >= days

    def append(self, from_table, to_table, pkey, **kwargs):
        '''
        make copy of input database table and import into another database table
        '''
        export_database, export_table = self.parse_database_table(from_table)
        import_database, import_table = self.parse_database_table(to_table)
        backup_directory = 'rethindb_export_' + export_database + '_' + export_table
        json_file = backup_directory + '/' + export_database + '/' + export_table + '.json'
        self.export_json(backup_directory, export_database, export_table, **kwargs)
        self.import_json(backup_directory, json_file, import_database, import_table, pkey, **kwargs)

    def parse_database_table(self, database_table):
        '''
        Parse the name of the database and table, return as tuple (database, table)
        '''
        database_table_split = database_table.split('.')
        if len(database_table_split) != 2:
            raise Exception('Name of database and table should be formatted as \'database.table\'')
        return database_table_split[0].strip().lower(), database_table_split[1].strip().lower()

    def export_json(self, backup_directory, export_database, export_table, rethink_host, auth_key, **kwargs):
        '''
        export export_database.export_table to backup_file
        '''
        print("Exporting database: " + export_database + ", table: " + export_table + ", to directory: " + backup_directory)
        if os.path.isdir(backup_directory):
            shutil.rmtree(backup_directory)
        if rethink_host!='localhost' and auth_key is not None:
            command = ['rethinkdb', 'export', '-c', rethink_host, '-a', auth_key, '-e', export_database + '.' + export_table, '-d', backup_directory, '--format', 'json']
        else:
            command = ['rethinkdb', 'export', '-e', export_database + '.' + export_table, '-d', backup_directory, '--format', 'json']
        try:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
        except:
            raise Exception("Couldn't export " + export_database + '.' + export_table + " to " + backup_directory)

    def import_json(self, backup_directory, json_file, import_database, import_table, pkey, rethink_host, auth_key, **kwargs):
        '''
        import the exported table into self.import_database from self.backup_file
        '''
        print("Importing into database: " + import_database + ", table: " + import_table + ", from file: " + json_file)
        if rethink_host!='localhost' and auth_key is not None:
            command = ['rethinkdb', 'import', '-c', rethink_host, '-a', auth_key, '--table', import_database+'.'+import_table, '-f', json_file, '--format', 'json', '--pkey', pkey, '--force']
        else:
            command = ['rethinkdb', 'import', '--table', import_database+'.'+import_table, '-f', json_file, '--format', 'json', '--pkey', pkey, '--force']
        try:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
            shutil.rmtree(backup_directory)
        except:
            raise Exception("Couldn't import into " + import_database + '.' + import_table + " from " + backup_directory+'/'+json_file)

    def sync_from_local(self, rethink_host, auth_key, export_database, export_table, pkey, **kwargs):
        '''
        Sync the local rethinkdb instance to an external rethinkdb instance
        Export documents in local database to external database
        '''
        print("Syncing local " + export_database+'.'+export_table + " table to external database")
        backup_directory = 'sync_from_local'
        json_file = backup_directory + '/' + export_database + '/' + export_table + '.json'
        print("Exporting documents from the local rethinkdb instance")
        self.export_json('localhost', None, backup_directory, export_database, export_table, **kwargs)
        print("Importing documents into the external rethinkdb instance")
        self.import_json(rethink_host, auth_key, backup_directory, json_file, pkey=pkey, **kwargs)

    def sync_to_local(self, rethink_host, auth_key, export_database, export_table, pkey, **kwargs):
        '''
        Sync the local rethinkdb instance to an external rethinkdb instance
        Export documents in external database to local database
        '''
        print("Syncing external " + export_database+'.'+export_table + " table to local database")
        backup_directory = 'sync_to_local'
        json_file = backup_directory + '/' + export_database + '/' + export_table + '.json'
        print("Exporting documents from the external rethinkdb instance")
        self.export_json(rethink_host, auth_key, backup_directory, export_database, export_table, **kwargs)
        print("Importing documents into the local rethinkdb instance")
        self.import_json('localhost', None, backup_directory, json_file, pkey=pkey, **kwargs)