import os, json
import rethinkdb as r
from upload import upload
from upload import get_parser

parser = get_parser()
parser.add_argument('--update_citations', default=False, action="store_true", help="update citation fields")
parser.add_argument('--update_locations', default=False, action="store_true", help="update location fields")
parser.add_argument('--update_passage_categories', default=False, action="store_true", help="update passage_category fields")
parser.add_argument('--update_groupings', default=False, action="store_true", help="update grouping fields")
parser.add_argument('--n_entrez', default=2500, type=int, help="number of entrez accessions to fetch at a time (default: 2500)")
parser.add_argument('--update_citations_tsv', default=None, type=str, help="update citations from tsv not entrez")


class update(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.location_fields = ['location', 'division', 'country', 'region']

    def update(self, update_citations, update_citations_tsv, update_locations, update_passage_categories, update_groupings, **kwargs):
        self.connect(**kwargs)
        if update_citations:
            self.update_citations(table=self.sequences_table, **kwargs)
        elif update_citations_tsv:
            self.update_citations_tsv(update_citations_tsv, table=self.sequences_table, **kwargs)
        elif update_locations:
            self.update_locations(table=self.viruses_table, **kwargs)
        elif update_passage_categories:
            self.update_passage_categories(table=self.sequences_table, **kwargs)
        elif update_groupings:
            self.update_groupings(self.viruses_table, self.sequences_table, **kwargs)
        else:
            self.update_genbank_documents(**kwargs)

    def update_citations(self, database, table, preview, index='accession', **kwargs):
        print("Updating citation fields")
        _, sequences = self.get_genbank_sequences(**kwargs)
        # the sequences come from the DB so we don't need to reformat them!
        # self.format_sequences(sequences, **kwargs)
        # self.match_duplicate_accessions(sequences, **kwargs)
        # self.match_database_duplicate_accessions(sequences, virus=self.virus, database=database)
        citation_keys = ['authors', 'title', 'journal', 'puburl', 'url', index]
        sequences = [{key: doc[key] for key in citation_keys} for doc in sequences]
        if not preview:
            print("Updating " + str(len(sequences)) + " sequence citations in " + database + "." + table)
            self.upload_to_rethinkdb(database, table, sequences, overwrite=True, index='accession')
        else:
            print("Preview of updates to be made, remove --preview to make updates to database")

    def get_genbank_sequences(self, email, **kwargs):
        if self.accessions is None:
            table = self.virus + "_sequences"
            accessions = self.get_accessions(self.database, table)
        else:
            accessions = [acc.strip() for acc in self.accessions.split(",")]
        self.entrez_email(email)
        gi = self.get_GIs(accessions, kwargs["n_entrez"])
        return self.get_entrez_viruses(gi, **kwargs)

    def get_accessions(self, database, table):
        '''
        Return documents from the table.database for which accession numbers are from genbank
        '''
        print("Getting accession numbers for sequences obtained from Genbank")
        accessions = list(r.db(database).table(table).filter((r.row["source"] == 'genbank') | (r.row["source"] == 'vipr')).get_field('accession').run())
        return accessions

    def update_locations(self, database, table, preview, **kwargs):
        print("Updating location fields")
        viruses = list(r.table(table).run())
        self.define_countries("source-data/geo_synonyms.tsv")
        self.define_regions("source-data/geo_regions.tsv")
        viruses = self.reassign_new_locations(viruses, self.location_fields, **kwargs)
        if not preview:
            print("Updating " + str(len(viruses)) + " virus locations in " + self.database + "." + self.viruses_table)
            del kwargs['overwrite']
            self.upload_to_rethinkdb(database, table, viruses, overwrite=True, index='strain')
        else:
            print("Preview of updates to be made, remove --preview to make updates to database")

    def reassign_new_locations(self, documents, location_fields, **kwargs):
        '''
        Use location fields to find location information that needs to be updated
        Uses the first location field that returns location information from self.determine_location
        Update location and region
        '''
        updated_documents = []
        fix_region = 'region' in location_fields # Do we also want to update region designations?
        if fix_region:
            location_fields.remove('region') # 'region' not a valid field for determine_location

        for doc in documents:
            old_location = tuple(doc[field] for field in location_fields if field in doc)
            for field in location_fields:
                if field in doc:
                    result = self.determine_location(self.snakecase_to_camelcase(doc[field]))
                    if result is not None:
                        if result != old_location:
                            doc['location'], doc['division'], doc['country'] = result
                            self.format_place(doc)
                            self.format_region(doc)
                            updated_documents.append(doc)
                    else:
                        print("couldn't parse %s for "%field, doc['strain'], self.snakecase_to_camelcase(doc[field]))
                    break  # found location information at detailed level

            if fix_region and doc not in updated_documents:
                try:
                    old_region = doc['region']
                    self.format_region(doc)
                    if doc['region'] != old_region:
                        updated_documents.append(doc)
                except Error as e:
                    print e
                    continue
        return updated_documents

    def update_passage_categories(self, **kwargs):
        print("No default passage category method defined, please use the specific virus update script, ie flu_update.py")

    def update_groupings(self, **kwargs):
        print("No default grouping method defined, please use the specific virus update script, ie flu_update.py")

    def update_genbank_documents(self, preview, **kwargs):
        '''
        Update all fields using genbank files
        '''
        print("Updating all fields from genbank")
        viruses, sequences = self.get_genbank_sequences(**kwargs)
        self.format_viruses(viruses, **kwargs)
        self.format_sequences(sequences, exclude_virus_methods=True, **kwargs)
        self.match_duplicate_strains(viruses, sequences, **kwargs)
        self.match_database_duplicate_strains(viruses, sequences, virus=self.virus, database=self.database, **kwargs)
        self.match_duplicate_accessions(sequences, **kwargs)
        self.match_database_duplicate_accessions(sequences, virus=self.virus, database=self.database)
        self.link_viruses_to_sequences(viruses, sequences)
        self.transfer_fields(viruses, sequences, self.virus_to_sequence_transfer_fields)
        if not preview:
            print("Updating viruses in " + self.database + "." + self.viruses_table)
            self.upload_to_rethinkdb(self.database, self.viruses_table, viruses, overwrite=True, index='strain')
            print("Updating sequences in " + self.database + "." + self.sequences_table)
            self.upload_to_rethinkdb(self.database, self.sequences_table, sequences, overwrite=True, index='accession')
        else:
            print("Viruses:")
            print(json.dumps(viruses[0], indent=1))
            print("Sequences:")
            print(json.dumps(sequences[0], indent=1))
            print("Preview of updates to be made, remove --preview to make updates to database")

    def update_citations_tsv(self, update_citations_tsv, database, table, preview, index='accession', **kwargs):
        print("Updating citation fields (manual table used for mumps)")
        citations = {}
        try:
            with open(update_citations_tsv, 'rU') as fh:
                for line in fh:
                    if not line.startswith('#'):
                        fields = line.strip().split("\t")
                        citations[fields[0]] = {
                            "authors": fields[1],
                            "title": fields[2],
                            "journal": fields[3],
                            "puburl": fields[4],
                            "url": fields[5]
                        }
        except IOError:
            raise Exception(update_citations_tsv, "not found")
        except IndexError:
            raise Exception(update_citations_tsv, "incorrect number of fields")

        sequences = list(r.db(database).table(table).run())
        sequences = [s for s in sequences if s["accession"] in citations.keys()]
        for s in sequences:
            if s["accession"] in citations:
                print("updating {} -> {}".format(s["accession"], citations[s["accession"]]["authors"]))
                for n in ["authors", "title", "journal", "puburl", "url"]:
                    s[n] = citations[s["accession"]][n]
            else:
                print("Sequence {} not in citations TSV".format(s["accession"]))

        if not preview:
            print("Updating " + str(len(sequences)) + " sequence citations in " + database + "." + table)
            self.upload_to_rethinkdb(database, table, sequences, overwrite=True, index='accession')
        else:
            print("Preview doesn't work")

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = update(**args.__dict__)
    connVDB.update(**args.__dict__)
