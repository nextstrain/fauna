import os, json, datetime
import rethinkdb as r
from Bio import SeqIO

class vdb_download(object):

    def __init__(self, database, virus, path='data/'):
        self.database = database
        self.virus_type = virus
        self.output_file_type = 'json'

        self.path = path
        self.current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))
        self.output_json_fname = self.path + self.virus_type + "/" + self.virus_type + "_" + self.current_date+ '.json'
        self.output_fasta_fname = self.path + self.virus_type + "/" + self.virus_type + "_" + self.current_date+ '.fasta'

        # connect to database
        try:
            r.connect(host="ec2-52-90-204-136.compute-1.amazonaws.com", port=28015, db=self.database, auth_key="KeHiybPoX8BM6aAhfYqy").repl()
            print("Connected to the \"" + self.database + "\" database")
        except:
            print("Failed to connect to the database, " + self.database)
            raise Exception

    def download_all_documents(self):
        '''
        download all documents from table
        :return:
        '''
        cursor = list(r.db(self.database).table(self.virus_type).run())
        print(type(cursor))
        for doc in cursor:
            self.pick_best_sequence(doc)
        self.write_json(cursor, self.output_json_fname)
        self.write_fasta(cursor, self.output_fasta_fname)

    def pick_best_sequence(self, document):
        '''
        find the best sequence in the given document. Currently by longest sequence.
        Resulting document is with flatter dictionary structure
        '''
        list_sequences = document['sequences']
        if len(list_sequences) == 1:
            best_sequence_info = document['sequences'][0]
        else:
            longest_sequence_pos = 0
            longest_sequence_length = len(document['sequences'][0]['sequence'])
            current_pos = 0
            for sequence_info in document['sequences']:
                if len(sequence_info['sequence']) > longest_sequence_length or (len(sequence_info['sequence']) ==
                                                        longest_sequence_length and sequence_info['accession'] is None):
                    longest_sequence_length = len(sequence_info['sequence'])
                    longest_sequence_pos = current_pos
                current_pos += 1
            best_sequence_info = document['sequences'][longest_sequence_pos]

        # create flatter structure for virus info
        for atr in best_sequence_info.keys():
            document[atr] = best_sequence_info[atr]
        del document['sequences']


    def write_json(self, data, file_name, indent=1):
        '''
        writes as list of viruses (dictionaries)
        '''
        print(type(file_name))
        print(file_name)
        try:
            handle = open(file_name, 'w')
        except:
            print("Couldn't open output file")
            print(file_name)
            raise FileNotFoundError
        else:
            json.dump(data, handle, indent=indent)
            handle.close()

    def write_fasta(self, viruses, file_name):
        fasta_fields = ['strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'location', 'source', 'locus']
        try:
            handle = open(file_name, 'w')
        except IOError:
            pass
        else:
            for v in viruses:
                handle.write(">")
                for field in fasta_fields:
                    handle.write(str(v[field]) + "|")
                handle.write("\n")
                handle.write(v['sequence'] + "\n")
            handle.close()


if __name__=="__main__":

    database = 'vdb'
    virus = 'Zika'

    upload_run = vdb_download(database, virus)
    upload_run.download_all_documents()