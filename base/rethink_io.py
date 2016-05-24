import rethinkdb as r

class rethink_io(object):
    def __init__(self, **kwargs):
        pass
    def connect_rethink(self, database, table, auth_key, rethink_host='localhost', **kwargs):
        '''
        Connect to rethink database,
        Check for existing table, otherwise create it
        '''
        if auth_key is None:
            try:
                r.connect(host=rethink_host, port=28015, db=database).repl()
                print("Connected to the \"" + database + "\" database")
            except:
                raise Exception("Failed to connect to the database, " + database)
        else:
            try:
                r.connect(host=rethink_host, port=28015, db=database, auth_key=auth_key).repl()
                print("Connected to the \"" + database + "\" database")
            except:
                raise Exception("Failed to connect to the database, " +database)

        existing_tables = r.db(database).table_list().run()
        if table not in existing_tables:
            raise Exception("No table exists yet for " + table + " available are " + str(existing_tables))