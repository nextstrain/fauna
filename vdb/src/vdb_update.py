import rethinkdb as r
from vdb_upload import vdb_upload
from vdb_upload import parser
from vdb_parse import vdb_parse
import re

class vdb_update(vdb_upload):
    def __init__(self, **kwargs):
        vdb_upload.__init__(self, **kwargs)
        vdb_parse.__init__(self, **kwargs)
        self.updateable_citation_fields = ['authors', 'title', 'url']
        self.updateable_sequence_fields = ['sequence']

    def update(self, **kwargs):
        if self.accessions is None:
            self.accessions = self.get_accessions()
        gi = self.get_GIs(self.accessions)
        self.viruses = self.get_entrez_viruses(gi, **kwargs)
        self.format()
        self.update_documents(**kwargs)

    def get_accessions(self):
        print("Getting accession numbers for sequences obtained from Genbank")
        cursor = list(r.db(self.database).table(self.virus).run())
        accessions = []
        for doc in cursor:
            for seq in doc['sequences']:
                if seq['source'] == 'Genbank':
                    accessions.append(seq['accession'])
        return accessions

    def update_documents(self, **kwargs):
        print("Checking for updates to " + str(len(self.viruses)) + " sequences")
        self.relaxed_strains()
        for virus in self.viruses:
            relaxed_name = self.relax_name(virus['strain'])
            if relaxed_name in self.strains:
                self.strain_name = self.strains[relaxed_name]
            else:
                self.strain_name = virus['strain']
            document = r.table(self.virus).get(self.strain_name).run()
            # Retrieve virus from table to see if it already exists
            if document is not None:
                self.updated = False
                self.update_sequence_citation_field(virus, document, 'accession', self.updateable_sequence_fields, self.updateable_citation_fields, **kwargs)

    def create_citations_field(self):
        cursor = list(r.db(self.database).table(self.virus).run())
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
            r.table(self.virus).get(strain).update({'sequences': sequence_info}).run()
            r.table(self.virus).get(strain).update({'citations': citation_info}).run()

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = vdb_update(**args.__dict__)
    connVDB.update(**args.__dict__)
