import rethinkdb as r
import datetime, os

class rethink_io(object):
    def __init__(self, **kwargs):
        pass

    def assign_rethink(self, rethink_host, auth_key, local=False, **kwargs):
        '''
        Assign rethink_host, auth_key, return as tuple
        '''
        if local:
            rethink_host = 'localhost'
        elif rethink_host is not None:
            rethink_host=rethink_host
        else:
            try:
                rethink_host = os.environ['RETHINK_HOST']
            except:
                raise Exception("Missing rethink host")
        if local:
            auth_key = None
        elif auth_key is not None:
            auth_key = auth_key
        else:
            try:
                auth_key = os.environ['RETHINK_AUTH_KEY']
            except:
                raise Exception("Missing rethink auth_key")
        return rethink_host, auth_key

    def connect_rethink(self, db, rethink_host='localhost', auth_key=None, **kwargs):
        '''
        Connect to rethink database,
        '''
        if rethink_host == 'localhost':
            try:
                conn = r.connect(host=rethink_host, port=28015, db=db).repl()
                print("Connected to the \"" + db + "\" database")
                return conn
            except:
                raise Exception("Failed to connect to the database, " + db)
        else:
            try:
                conn = r.connect(host=rethink_host, port=28015, db=db, auth_key=auth_key).repl()
                print("Connected to the \"" + db + "\" database")
                return conn
            except:
                raise Exception("Failed to connect to the database, " + db)

    def check_table_exists(self, database, table):
        '''
        Check for existing table
        '''
        existing_tables = r.db(database).table_list().run()
        if table not in existing_tables:
            raise Exception("No table exists yet for " + table + " available are " + str(existing_tables))

    def get_upload_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y-%m-%d'))

    def get_upload_timestamp(self):
        return str(datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y-%m-%d-%H-%M'))

    def check_optional_attributes(self, doc, optional_fields):
        '''
        Reassign unknowns from '?' or '' to 'None'
        Create and assign 'None' to optional attributes that don't exist
        '''
        for key in doc.keys():
            if doc[key] == '?' or doc[key] == '':
                doc[key] = None
            if isinstance(doc[key], (str)):
                doc[key] = doc[key].strip()
        for atr in optional_fields:
            if atr not in doc:
                doc[atr] = None

    def check_required_attributes(self, doc, required_fields, index_fields):
        '''
        Checks that required upload attributes are present and not equal to None for given virus
        :return: returns true if it has all required upload attributes, else returns false and prints missing attributes
        '''
        missing_attributes = []
        for atr in required_fields:
            if atr in doc and doc[atr] is not None:
                pass
            else:
                missing_attributes.append(atr)
        if len(missing_attributes) > 0:
            print("This document is missing a required attribute and will be removed from upload sequences")
            print([doc[index] for index in index_fields if index in doc])
            print("Missing attributes: " + str(missing_attributes))
            return False
        else:
            return True

    def delete_extra_fields(self, doc, fields, index_fields):
        '''
        Delete attributes not defined as a field
        '''
        for key in doc.keys():
            if key not in fields:
                print("Deleting document info " + key + ": " + doc[key] + " from " + str([doc[index] for index in index_fields if index in doc]))
                del doc[key]