import rethinkdb as r
import datetime

class rethink_io(object):
    def __init__(self, **kwargs):
        pass

    def connect_rethink(self, database, table, auth_key, rethink_host='localhost', **kwargs):
        '''
        Connect to rethink database,
        Check for existing table
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


    def get_upload_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))

    def check_optional_attributes(self, documents, optional_fields):
        '''
        Reassign unknowns from '?' or '' to 'None'
        Create and assign 'None' to optional attributes that don't exist
        '''
        for doc in documents:
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