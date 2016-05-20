import os, argparse, shutil
import rethinkdb as r

parser = argparse.ArgumentParser()
#parser.add_argument('-db', '--database', default='vdb', help="database to make copy of and import to")
#parser.add_argument('-dv', '--virus', default='zika', help="virus to make copy of and import to")
parser.add_argument('-edb', '--export_database', default='vdb', help="database to make copy of")
parser.add_argument('-ev', '--export_virus', default='zika', help="virus table to copy")
parser.add_argument('-idb', '--import_database', default='test_vdb', help="database that will be changed")
parser.add_argument('-iv', '--import_virus', default='zika', help="virus table to import to")
parser.add_argument('--copy', default=False, action="store_true", help="copy database and table to another database and table")
parser.add_argument('--import_json', default=False, action="store_true", help="upload to import database and import virus from json")
parser.add_argument('--fname', default=None, help="json to import from")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")


class vdb_copy(object):
    def __init__(self, export_database, export_virus, import_database, import_virus, rethink_host=None, auth_key=None, **kwargs):
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

    def copy(self, export_database, export_virus, **kwargs):
        '''
        make copy of input database table and import into another database table
        '''
        backup_file = 'rethindb_export_' + export_database + '_' + export_virus
        self.export_table(backup_file, export_database, export_virus, **kwargs)
        self.import_table(backup_file, export_database, export_virus, **kwargs)

    def export_table(self, backup_file, export_database, export_virus, **kwargs):
        '''
        backup self.export_database self.virus table to self.backup_file
        '''
        print("Making backup of database: " + export_database.lower() + ", table: " + export_virus.lower() + ", to file: " + backup_file)
        command = 'rethinkdb export -c ' + self.rethink_host + ' -a ' + self.auth_key + ' -e ' + export_database.lower() + '.' + export_virus.lower() + ' -d ' + backup_file + ' --format json'
        os.system(command)

    def import_table(self, backup_file, export_database, export_virus, import_database, import_virus, **kwargs):
        '''
        import the exported table into self.import_database from self.backup_file
        '''
        print("Importing into database: " + import_database + ", table: " + import_virus + ", from file: " + backup_file)
        json_file = backup_file + '/' + export_database + '/' + export_virus + '.json'
        command = 'rethinkdb import -c ' + self.rethink_host + ' -a ' + self.auth_key + ' --table ' + import_database + '.' + import_virus + ' -f ' + json_file + ' --format json --pkey strain --force'
        os.system(command)
        shutil.rmtree(backup_file)

    def import_json(self, fname, import_database, import_virus, **kwargs):
        '''
        import the json file into self.import_database from fname
        '''
        print("Importing into database: " + import_database + ", table: " + import_virus + ", from file: " + fname)
        command = 'rethinkdb import -c ' + self.rethink_host + ' -a ' + self.auth_key + ' --table ' + import_database + '.' + import_virus + ' -f ' + fname + ' --format json --pkey strain --force'
        os.system(command)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = vdb_copy(**args.__dict__)
    if args.import_json:
        connVDB.import_json(**args.__dict__)
    else:
        connVDB.copy(**args.__dict__)

