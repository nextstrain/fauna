import os, argparse, sys, time, datetime
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from base.rethink_interact import rethink_interact
from vdb.restore import parser

class restore(object):
    def __init__(self, database, **kwargs):
        self.upload_hour = 3
        self.database = database.lower()
        if 'database' in kwargs:
            kwargs['database'] = self.database
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_interact = rethink_interact()

    def restore(self, virus, restore_date, **kwargs):
        '''
        restore vdb to previous version
        '''
        kwargs['rethink_host'] = self.rethink_host
        kwargs['auth_key'] = self.auth_key
        virus = virus.lower()
        self.rethink_interact.restore(restore_table=virus, restore_date=restore_date, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = restore(**args.__dict__)
    connVDB.restore(**args.__dict__)