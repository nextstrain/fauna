import rethinkdb as r
from upload import upload
from upload import parser

class update(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def update(self, email, **kwargs):
        if self.accessions is None:
            table = self.virus + "_sequences"
            accessions = self.get_accessions(self.database, table)
        else:
            accessions = [acc.strip() for acc in self.accessions.split(",")]
        self.entrez_email(email)
        gi = self.get_GIs(accessions)
        viruses, sequences = self.get_entrez_viruses(gi, **kwargs)
        self.format(viruses, **kwargs)
        self.format(sequences, exclude_virus_methods=True, **kwargs)
        self.link_viruses_to_sequences(viruses, sequences)
        self.transfer_fields(viruses, sequences, self.virus_to_sequence_transfer_fields)
        print("Updating viruses in " + self.database + "." + self.viruses_table)
        self.upload_documents(self.viruses_table, viruses, 'strain', **kwargs)
        print("Updating sequences in " + self.database + "." + self.sequences_table)
        self.upload_documents(self.sequences_table, sequences, 'accession', **kwargs)

    def get_accessions(self, database, table):
        print("Getting accession numbers for sequences obtained from Genbank")
        cursor = list(r.db(database).table(table).run())
        accessions = []
        for doc in cursor:
            if 'source' in doc:
                if doc['source'].lower() == 'genbank' or doc['source'].lower() == 'vipr':
                    accessions.append(doc['accession'])
        return accessions

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = update(**args.__dict__)
    connVDB.update(**args.__dict__)
