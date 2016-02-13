

class vdb_virus(object):

    def __init__(self, strain, date, country, list_sequences, division=None, location=None,meta=None):
        self.strain = strain
        self.date = date
        self.country = country
        self.sequences = list_sequences  # list of dictionaries containing sequence information
        # check that sequences contains (Accession, Database, Region, Sequence)
        self.division = division
        self.location = location
        self.meta = meta  # dictionary containing meta information

    def fasta_format(self):
        '''

        :return: string containing virus information in fasta format
        '''

    def json_format(self):
        '''

        :return: string containing virus information in json format
        '''