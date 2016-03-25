import os, argparse, shutil
import rethinkdb as r

parser = argparse.ArgumentParser()
parser.add_argument('-edb', '--export_database', default='vdb', help="database to make backup of")
parser.add_argument('-idb', '--import_database', default='test', help="database that will be restored to backup copy")
parser.add_argument('-v', '--virus', help="virus table to backup or copy")
parser.add_argument('--fname', help="input file name")
parser.add_argument('--path', default=None, help="path to fasta file, default is \"data/virus/\"")
parser.add_argument('--host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")


class vdb_copy(object):
    def __init__(self, **kwargs):
        if 'virus' in kwargs:
            self.virus = kwargs['virus'].title()
        if 'export_database' in kwargs:
            self.export_database = kwargs['export_database']
        if 'import_database' in kwargs:
            self.import_database = kwargs['import_database']
        self.backup_file = 'rethindb_export_' + self.export_database + '_' + self.virus


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
        print("Making backup of database: " + self.export_database + ", table: " + self.virus + ", to file: " + self.backup_file)
        command = 'rethinkdb export -c ' + self.host + ' -a ' + self.auth_key + ' -e ' + self.export_database + '.' + self.virus + ' -d ' + self.backup_file + ' --format json'
        os.system(command)

    def import_table(self):
        '''
        import the exported table into self.import_database from self.backup_file
        '''
        print("Restoring database: " + self.import_database + ", table: " + self.virus + ", from file: " + self.backup_file)
        json_file = self.backup_file + '/' + self.export_database + '/' + self.virus + '.json'
        command = 'rethinkdb import -c ' + self.host + ' -a ' + self.auth_key + ' --table ' + self.import_database + '.' + self.virus + ' -f ' + json_file + ' --format json --pkey strain --force'
        os.system(command)
        shutil.rmtree(self.backup_file)

if __name__=="__main__":
    args = parser.parse_args()
    run = vdb_copy(**args.__dict__)
    run.copy()


