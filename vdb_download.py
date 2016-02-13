import os
from pymongo import MongoClient


class vdb_download(object):

    def __init__(self, collection):
        client = MongoClient()
        self.db = client['vdb']  # connect to specified database, database must be string
        self.collection = self.db[collection]  # connect to specified collection, collection must be string

    def download_all_documents(self):
        '''
        download all documents from vdb
        :return:
        '''

    def output_result(self):
        '''
        Output to fasta or json document
        :return:
        '''