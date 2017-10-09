import os, json, sys
import random
import argparse
import rethinkdb as r
from parse import parse
sys.path.append('')
from base.rethink_io import rethink_io

parser = argparse.ArgumentParser()
parser.add_argument('--path', default='input_data/', help="Path to directory containing upload file, default is input_data/")
parser.add_argument('--preview', default=False, action="store_true",  help ="If included, preview a virus document to be uploaded")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('filename', help="file to upload to rethink database")

class json_upload(parse):

    def __init__(self, filename, **kwargs):
        parse.__init__(self, **kwargs)
        self.uploadable_databases = ['sequence', 'virus', 'titer']
        self.read_sacra_json(filename, **kwargs)

    def read_sacra_json(self, filename, **kwargs):
        '''
        Parse header material from a Sacra JSON to determine
        the upload loation and dataset for the upload run.
        '''
        with open(filename, 'r') as f:
            dataset = json.load(f)
            self.datatype = str(dataset["datatype"])
            self.virus = str(dataset["virus"])
            # self.virus_set = dataset["viruses"] TODO: Add the logical set of viruses to sacra build
            self.data = dataset["data"]

    def upload(self, preview=False, **kwargs):
        '''
        Upload self.data to database
        '''
        self.connect(**kwargs)
        print("Uploading Viruses to VDB")
        # viruses = self.virus_set
        if not preview:
            # TODO: Deal with viruses once Sacra is updated
            # print("Uploading viruses to virus." + self.virus)
            # self.upload_documents(self.viruses_table, self.virus_set, index='strain', **kwargs)
            print("Uploading sequences to " + self.datatype + "." + self.virus)
            self.upload_documents(self.virus, self.data.values(), self.datatype, index='accession', **kwargs)
        else:
            # print("Viruses:")
            # print(json.dumps(viruses[0], indent=1))
            print("Sequences:")
            print(json.dumps(random.choice(list(self.data.items())), indent=1))
            print("Remove \"--preview\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def upload_documents(self, table, documents, database, replace=False, **kwargs):
        if replace:
            print("Deleting documents in database:" + database + "." + table)
            r.table(table).delete().run()
        print("Inserting ", len(documents), "documents")
        self.upload_to_rethinkdb(database, table, documents, **kwargs)

    def upload_to_rethinkdb(self, database, table, documents, overwrite=False, optimal_upload=200, **kwargs):
        if len(documents) > optimal_upload:
            list_documents = [documents[x:x+optimal_upload] for x in range(0, len(documents), optimal_upload)]
        else:
            list_documents = [documents]
        print("Uploading to rethinkdb in " + str(len(list_documents)) + " batches of " + str(optimal_upload) + " documents at a time")
        inserted = 0
        replaced = 0
        for list_docs in list_documents:
            if not overwrite:
                document_changes = r.table(table).insert(list_docs, conflict=lambda id, old_doc, new_doc: rethinkdb_updater(id, old_doc, new_doc), return_changes=True).run()
            else:
                document_changes = r.table(table).insert(list_docs, conflict=lambda id, old_doc, new_doc: rethinkdb_updater_overwrite(id, old_doc, new_doc), return_changes=True).run()
            self.update_timestamp(table, document_changes, **kwargs)
            if document_changes['errors']>0:
                print("Errors were made when inserting the documents", document_changes['errors'])
                print(document_changes['errors'])
            inserted += document_changes['inserted']
            replaced += document_changes['replaced']
        print("Ended up inserting " + str(inserted) + " documents into " + database + "." + table)
        print("Ended up updating " + str(replaced) + " documents in " + database + "." + table)

    def connect(self, virus_upload=False, **kwargs):
        if self.datatype not in self.uploadable_databases:
            raise Exception("Can't upload to this database: " + self.datatype, "add to list of databases allowed", self.uploadable_databases)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        if not virus_upload:
            self.rethink_io.connect_rethink(self.datatype, self.rethink_host, self.auth_key)
            self.rethink_io.check_table_exists(self.datatype, self.virus)
        else:
            self.rethink_io.connect_rethink('virus', self.rethink_host, self.auth_key)
            self.rethink_io.check_table_exists('virus', self.virus)

    def update_timestamp(self, table, document_changes, index, **kwargs):
        '''
        Update the timestamp field in the rethink table if changes have been made to the documents
        '''
        updated_documents = []
        if 'changes' in document_changes:
            for doc in document_changes['changes']:
                if doc['new_val'] is not None:
                    updated_documents.append({index: doc['new_val'][index], 'timestamp': self.rethink_io.get_upload_timestamp()})
        if len(updated_documents) > 0:
            r.table(table).insert(updated_documents, conflict='update').run()

##### End class ######

def rethinkdb_updater(id, old_doc, new_doc):
    return (new_doc.keys().set_union(old_doc.keys()).map(lambda key:
        r.branch(old_doc.keys().contains(key).and_(new_doc.keys().contains(key).not_()),
            [key, old_doc[key]],
            new_doc.keys().contains(key).and_(old_doc.keys().contains(key).not_()),
            [key, new_doc[key]],
            r.branch(key.eq('sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences'])],
                key.eq('number_sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences']).count()],
                key.eq('timestamp').or_(key.eq('virus_inclusion_date')).or_(key.eq('sequence_inclusion_date')),
                [key, old_doc[key]],
                r.branch(old_doc[key].eq(None).and_(new_doc[key].eq(None).not_()),
                    [key, new_doc[key]],
                    [key, old_doc[key]])
            )
        )
    )).coerce_to('object')

def rethinkdb_updater_overwrite(id, old_doc, new_doc):
    return (new_doc.keys().set_union(old_doc.keys()).map(lambda key:
        r.branch(old_doc.keys().contains(key).and_(new_doc.keys().contains(key).not_()),
            [key, old_doc[key]],
            new_doc.keys().contains(key).and_(old_doc.keys().contains(key).not_()),
            [key, new_doc[key]],
            r.branch(key.eq('sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences'])],
                key.eq('number_sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences']).count()],
                key.eq('timestamp').or_(key.eq('virus_inclusion_date')).or_(key.eq('sequence_inclusion_date')),
                [key, old_doc[key]],
                [key, new_doc[key]]
            )
        )
    )).coerce_to('object')

if __name__=="__main__":
    args = parser.parse_args()
    J = json_upload(**args.__dict__)
    J.upload(**args.__dict__)
