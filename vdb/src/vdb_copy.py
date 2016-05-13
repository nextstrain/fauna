import os, argparse, shutil
import rethinkdb as r

parser = argparse.ArgumentParser()
parser.add_argument('-edb', '--export_database', default='vdb', help="database to make copy of")
parser.add_argument('-ev', '--export_virus', default='zika', help="virus table to copy")
parser.add_argument('-idb', '--import_database', default='test_vdb', help="database that will be changed")
parser.add_argument('-iv', '--import_virus', default='zika', help="virus table to import to")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")


class vdb_copy(object):
    def __init__(self, export_database, export_virus, import_database, import_virus, rethink_host=None, auth_key=None, **kwargs):
        self.export_virus = export_virus.lower()
        self.export_database = export_database.lower()
        self.import_database = import_database.lower()
        self.import_virus = import_virus.lower()
        self.backup_file = 'rethindb_export_' + self.export_database + '_' + self.export_virus
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

    def copy(self):
        '''
        make copy of input database table and import into another database table
        '''
        self.export_table()
        self.import_table()

    def export_table(self):
        '''
        backup self.export_database self.virus table to self.backup_file
        '''
        print("Making backup of database: " + self.export_database + ", table: " + self.export_virus + ", to file: " + self.backup_file)
        command = 'rethinkdb export -c ' + self.rethink_host + ' -a ' + self.auth_key + ' -e ' + self.export_database + '.' + self.export_virus + ' -d ' + self.backup_file + ' --format json'
        os.system(command)

    def import_table(self):
        '''
        import the exported table into self.import_database from self.backup_file
        '''
        print("Importing into database: " + self.import_database + ", table: " + self.import_virus + ", from file: " + self.backup_file)
        json_file = self.backup_file + '/' + self.export_database + '/' + self.export_virus + '.json'
        command = 'rethinkdb import -c ' + self.rethink_host + ' -a ' + self.auth_key + ' --table ' + self.import_database + '.' + self.import_virus + ' -f ' + json_file + ' --format json --pkey strain --force'
        os.system(command)
        shutil.rmtree(self.backup_file)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = vdb_copy(**args.__dict__)
    connVDB.copy()


