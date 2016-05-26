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
            command = ['rethinkdb', 'dump', '-c', rethink_host, '-e', database + '.' + dump_table, '-f', dump_file]
        try:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(command, stdout=devnull, stderr=subprocess.STDOUT)
        except:
            raise Exception("Couldn't dump tar file, make sure " + dump_file + " doesn't exist")

    def restore(self, database, restore_table, restore_date, rethink_host='localhost', auth_key=None, **kwargs):
        '''
        '''
        if restore_table is not None and restore_date is not None:
            print("Restoring database " + database + " and table " + restore_table.lower() + " to date " + restore_date)
            restore_location = database+'.'+restore_table
            restore_file = restore_date + "_" + database + "_" + restore_table.lower() + '.tar.gz'
            self.get_file(fname=restore_file, **kwargs)
            if rethink_host != 'localhost' and auth_key is not None:
                command = ['rethinkdb', 'restore', restore_file, '-i', restore_location, '-c', rethink_host, '-a', auth_key, '--force']
            else:
                command = ['rethinkdb', 'restore', restore_file, '-i', restore_location, '-c', rethink_host, '--force']
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