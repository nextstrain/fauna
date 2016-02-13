import os
from pymongo import MongoClient


class vdb_upload(object):

    def __init__(self, collection):
        client = MongoClient()
        self.db = client['vdb']  # connect to specified database, database must be string
        self.collection = self.db[collection]  # connect to specified collection, collection must be string



    def parse_fasta(self):
        '''
        parse input fasta file, create list of documents
        :return: list of documents(dictionaries of attributes) to upload
        '''


    def upload_documents(self, documents):
        '''

        :return:
        '''

        self.collection.insert_many