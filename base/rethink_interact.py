import os, shutil, sys, datetime, re, subprocess, json
from rethinkdb import r
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

    def backup_s3(self, database, path='', delete_expired=False, **kwargs):
        '''
        make backup of every table in database, upload to s3 bucket
        '''
        print("Backing up " + database + " on " + self.rethink_io.get_upload_date())
        if not os.path.isdir(path):
            os.makedirs(path)
        tables = r.db(database).table_list().run()
        bucket = self.connect_S3(**kwargs)
        for table in tables:
            dump_file = self.rethink_io.get_upload_date() + '_' + database + '_' + table + '.tar.gz'
            self.dump(database=database, dump_table=table, dump_file=dump_file, **kwargs)
            shutil.move(dump_file, path+'/'+dump_file)
            bucket.upload_file(path+'/'+dump_file, dump_file)
        print("Successfully backed up")
        if delete_expired:
            self.delete_expired_s3_backups(bucket, **kwargs)

    def backup_local(self, database, path='', delete_expired=False, **kwargs):
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
        if delete_expired:
            self.delete_expired_local_backups(path=path, **kwargs)

    def dump(self, database, dump_table, dump_file, rethink_host='localhost', auth_key=None, **kwargs):
        '''
        backup self.export_database table to self.backup_file
        :return:
        '''
        if rethink_host!='localhost':
            command = ['rethinkdb', 'dump', '-c', rethink_host, '-a', auth_key, '-e', database + '.' + dump_table, '-f', dump_file]
        else:
            command = ['rethinkdb', 'dump', '-e', database + '.' + dump_table, '-f', dump_file]
        try:
            with open(os.devnull, 'wb') as devnull:
                print " ".join(command)
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT, shell=True)
        except:
            raise Exception("Couldn't dump tar file, make sure " + dump_file + " doesn't exist")

    def restore(self, database, restore_table, restore_date, rethink_host='localhost', auth_key=None, **kwargs):
        '''
        Database is database to restore to
        Backup_database is database from which the backup file was created
        '''
        if restore_table is not None and restore_date is not None:
            restore_file = restore_date + "_" + database + "_" + restore_table + '.tar.gz'
            restore_to_location = database+'.'+restore_table
            print("Restoring database " + database + "." + restore_table + " to date " + restore_date + " from file " + restore_file)
            self.get_file(fname=restore_file, **kwargs)
            if rethink_host != 'localhost' and auth_key is not None:
                command = ['rethinkdb', 'restore', restore_file, '-i', restore_to_location, '-c', rethink_host, '-a', auth_key, '--force']
            else:
                command = ['rethinkdb', 'restore', restore_file, '-i', restore_to_location, '--force']
            try:
                with open(os.devnull, 'wb') as devnull:
                    subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
                os.remove(restore_file)
            except:
                raise Exception("Couldn't restore " + str(restore_to_location) + " with file " + str(restore_file))
        else:
            raise Exception("define arguments \'--restore_table\' and \'--restore_date\'")

    def get_file(self, backup_s3, backup_local, fname, path='', **kwargs):
        '''
        Download file from s3 to current directory
        '''
        if backup_s3:
            bucket = self.connect_S3(**kwargs)
            existing_files = [str(object.key) for object in bucket.objects.all()]
            if fname in existing_files:
                bucket.download_file(fname, fname)
            else:
                raise Exception("Specified file could not be found in s3 bucket", fname)
        elif backup_local:
            shutil.copy(path+'/'+fname, fname)
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

    def delete_expired_local_backups(self, path, days_to_expiration, **kwargs):
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

    def append(self, from_table, to_table, rethink_host, auth_key, **kwargs):
        '''
        make copy of input database table and import into another database table
        '''
        export_database, export_table = self.parse_database_table(from_table)
        import_database, import_table = self.parse_database_table(to_table)
        directory = 'append'
        json_file = directory + '/' + export_database + '_' + export_table + '.json'

        conn = self.rethink_io.connect_rethink(export_database, rethink_host, auth_key)
        self.rethink_io.check_table_exists(export_database, export_table)
        self.export_json(directory, json_file, export_database, export_table, **kwargs)
        conn.close()

        conn = self.rethink_io.connect_rethink(import_database, rethink_host, auth_key)
        self.rethink_io.check_table_exists(import_database, import_table)
        self.import_json(directory, json_file, import_database, import_table, **kwargs)
        conn.close()

    def parse_database_table(self, database_table):
        '''
        Parse the name of the database and table, return as tuple (database, table)
        '''
        database_table_split = database_table.split('.')
        if len(database_table_split) != 2:
            raise Exception('Name of database and table should be formatted as \'database.table\'')
        return database_table_split[0].strip().lower(), database_table_split[1].strip().lower()

    def export_json(self, directory, json_file, export_database, export_table, **kwargs):
        '''
        export export_database.export_table to backup_file
        '''
        print("Exporting database: " + export_database + ", table: " + export_table + ", to file: " + json_file)
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
        try:
            documents = list(r.db(export_database).table(export_table).run())
            write_json(documents, json_file)
        except:
            raise Exception("Couldn't export " + export_database + '.' + export_table + " to " + directory)

    def import_json(self, directory, json_file, import_database, import_table, **kwargs):
        '''
        import the exported table into self.import_database from self.backup_file
        '''
        print("Importing into database: " + import_database + ", table: " + import_table + ", from file: " + json_file)
        try:
            documents = read_json(json_file)
            self.sync_via_timestamp(import_table, documents, **kwargs)
            shutil.rmtree(directory)
        except:
            raise Exception("Couldn't import into " + import_database + '.' + import_table + " from " + json_file)

    def sync_via_timestamp(self, table, documents, key='strain', **kwargs):
        '''
        '''
        print(key)
        for document in documents:
            result = r.table(table).get(document[key]).run()
            if result is None:
                r.table(table).insert(document).run()
            else:
                if document['timestamp'] > result['timestamp']:
                    r.table(table).insert(document, conflict='replace').run()

    def push(self, local_table, remote_table, rethink_host, auth_key, **kwargs):
        '''
        push local database documents to remote database
        '''
        export_database, export_table = self.parse_database_table(local_table)
        import_database, import_table = self.parse_database_table(remote_table)
        print("Syncing local " + export_database+'.'+export_table + " table to external database")
        directory = 'push'
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
        json_file = directory + '/' + export_database + '_' + export_table + '.json'

        print("Exporting documents from the local rethinkdb instance")
        conn = self.rethink_io.connect_rethink(export_database, rethink_host='localhost', auth_key=None)
        self.rethink_io.check_table_exists(export_database, export_table)
        self.export_json(directory, json_file, export_database, export_table)
        conn.close()

        print("Importing documents into the external rethinkdb instance")
        self.rethink_io.connect_rethink(import_database, rethink_host, auth_key)
        self.rethink_io.check_table_exists(import_database, import_table)
        self.import_json(directory, json_file, import_database, import_table)

    def pull(self, local_table, remote_table, rethink_host, auth_key, **kwargs):
        '''
        pull remote database documents to local database
        '''
        export_database, export_table = self.parse_database_table(remote_table)
        import_database, import_table = self.parse_database_table(local_table)
        print("Syncing external " + export_database+'.'+export_table + " table to local database")
        directory = 'pull'
        if os.path.isdir(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)
        json_file = directory + '/' + export_database + '_' + export_table + '.json'

        print("Exporting documents from the external rethinkdb instance")
        conn = self.rethink_io.connect_rethink(export_database, rethink_host, auth_key)
        self.rethink_io.check_table_exists(export_database, export_table)
        self.export_json(directory, json_file, export_database, export_table)
        conn.close()

        print("Importing documents into the local rethinkdb instance")
        self.rethink_io.connect_rethink(import_database, rethink_host='localhost', auth_key=None)
        self.rethink_io.check_table_exists(import_database, import_table)
        self.import_json(directory, json_file, import_database, import_table)


def read_json(file_name):
    try:
        handle = open(file_name, 'r')
        data = json.load(handle)
        handle.close()
        return data
    except:
        raise Exception("Couldn't read json " + file_name)

def write_json(data, file_name, indent=1):
    try:
        handle = open(file_name, 'w')
        json.dump(data, handle, indent=indent)
        handle.write("\n")
        handle.close()
    except:
        raise Exception("Couldn't write to json " + file_name)
