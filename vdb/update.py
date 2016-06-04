import rethinkdb as r
from upload import upload
from upload import parser

class update(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.updateable_citation_fields = ['authors', 'title', 'url']
        self.updateable_sequence_fields = ['sequence']

    def update(self, **kwargs):
        if self.accessions is None:
            accessions = self.get_accessions()
        else:
            accessions = [acc.strip() for acc in self.accessions.split(",")]
        gi = self.get_GIs(accessions)
        self.viruses = self.get_entrez_viruses(gi, **kwargs)
        self.format()
        self.format_schema()
        self.update_documents(**kwargs)

    def get_accessions(self):
        print("Getting accession numbers for sequences obtained from Genbank")
        cursor = list(r.db(self.database).table(self.table).run())
        accessions = []
        for doc in cursor:
            index = 0
            for seq in doc['citations']:
                if 'source' in seq:
                    if seq['source'].lower() == 'genbank' or seq['source'].lower() == 'vipr':
                        accessions.append(doc['sequences'][index]['accession'])
                index += 1
        return accessions

    def update_documents(self, **kwargs):
        print("Checking for updates to " + str(len(self.viruses)) + " sequences")
        self.relaxed_strains()
        db_relaxed_strains = self.relaxed_strains()
        for virus in self.viruses:
            relaxed_name = virus['strain']
            if self.relax_name(virus['strain']) in db_relaxed_strains:
                relaxed_name = db_relaxed_strains[self.relax_name(virus['strain'])]
            try:
                document = r.table(self.table).get(relaxed_name).run()
            except:
                print(virus)
                raise Exception("Couldn't retrieve this virus")
            # Retrieve virus from table to see if it already exists
            if document is not None:
                updated = self.update_sequence_citation_field(document, virus, 'accession', self.updateable_sequence_fields, self.updateable_citation_fields, **kwargs)
                if updated:
                    document['timestamp'] = virus['timestamp']
                    r.table(self.table).insert(document, conflict="replace").run()

    def create_citations_field(self):
        cursor = list(r.db(self.database).table(self.table).run())
        for doc in cursor:
            strain = doc['strain']
            sequence_fields = self.sequence_upload_fields + self.sequence_optional_fields
            sequence_info = [{}]
            citation_info = [{}]
            for atr in doc['sequences'][0]:
                if atr not in sequence_fields:
                    citation_info[0][atr] = doc['sequences'][0][atr]
                else:
                    sequence_info[0][atr] = doc['sequences'][0][atr]
            for field in self.citation_optional_fields:
                if field not in citation_info[0]:
                    citation_info[0][field] = '?'
            r.table(self.table).get(strain).update({'sequences': sequence_info}).run()
            r.table(self.table).get(strain).update({'citations': citation_info}).run()

    def add_attribute(self):
        cursor = list(r.db(self.database).table(self.table).run())
        for doc in cursor:
            strain = doc['strain']
            r.table(self.table).get(strain).update({'host': 'human'}).run()
            r.table(self.table).get(strain).update({'subtype': None}).run()

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = update(**args.__dict__)
    connVDB.update(**args.__dict__)
