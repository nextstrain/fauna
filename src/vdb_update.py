import rethinkdb as r
from vdb_upload import vdb_upload
from vdb_upload import parser
from vdb_parse import vdb_parse

class vdb_update(vdb_upload):
    def __init__(self, **kwargs):
        vdb_upload.__init__(self, **kwargs)
        vdb_parse.__init__(self, **kwargs)
        self.updateable_fields = ['authors', 'title', 'url', 'sequence']

    def update(self):
        accessions = self.get_accessions()
        self.viruses = self.access_entrez(accessions)
        self.format()
        self.update_documents()

    def get_accessions(self):
        print("Getting accession numbers for sequences from Genbank")
        cursor = list(r.db(self.database).table(self.virus).run())
        accessions = []
        for doc in cursor:
            for seq in doc['sequences']:
                if seq['source'] == 'Genbank':
                    accessions.append(seq['accession'])
        return accessions

    def update_documents(self):
        print("Checking for updates to " + str(len(self.viruses)) + " sequences")
        for virus in self.viruses:
            # Retrieve virus from table to see if it already exists
            document = r.table(self.virus).get(virus['strain']).run()
            if document is not None:
                self.updated = False
                self.update_sequence_field(virus, document, 'accession')

    def update_sequence_field(self, virus, document, check_field):
        '''
        Checks for matching viruses by comparing sequence and accession, updates other attributes if needed
        '''
        doc_seqs = document['sequences']
        virus_seq = virus['sequences'][0]
        updated_sequence = False
        for doc_sequence_info in doc_seqs:
            if doc_sequence_info[check_field] == virus_seq[check_field]:
                for field in self.updateable_fields:
                    if field not in doc_sequence_info:
                        doc_sequence_info[field] = virus_seq[field]
                        print("Creating sequences field ", field, " assigned to ", virus_seq[field])
                        updated_sequence = True
                    elif (field in virus_seq and virus_seq[field] is not None and doc_sequence_info[field] != virus_seq[field]):
                        print("Updating virus field " + str(field) + ", from " + str(doc_sequence_info[field]) + " to " + str(virus_seq[field]))
                        doc_sequence_info[field] = virus_seq[field]
                        updated_sequence = True
        if updated_sequence:
            r.table(self.virus).get(virus['strain']).update({"sequences": doc_seqs}).run()
            r.table(self.virus).get(virus['strain']).update({'date_modified': virus['date_modified']}).run()
            self.updated = True

if __name__=="__main__":
    args = parser.parse_args()
    run = vdb_update(**args.__dict__)
    run.update()