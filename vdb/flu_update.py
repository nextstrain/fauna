import os
import numpy as np
from rethinkdb import r
from Bio import SeqIO
from update import update
from flu_upload import flu_upload
from update import parser

class flu_update(update, flu_upload):
    def __init__(self, **kwargs):
        update.__init__(self, **kwargs)
        flu_upload.__init__(self, **kwargs)
        self.outgroups = {lineage: SeqIO.read('source-data/'+lineage+'_outgroup.gb', 'genbank') for lineage in ['H3N2', 'H1N1', 'H1N1pdm', 'Vic', 'Yam']}
        self.outgroup_patterns = {'H3N2': ('a', 'h3n2', 'seasonal_h3n2'),
                                  'H1N1': ('a', 'h1n1', 'seasonal_h1n1'),
                                  'H1N1pdm': ('a', 'h1n1', 'seasonal_h1n1pdm'),
                                  'Vic': ('b', None, 'seasonal_vic'),
                                  'Yam': ('b', None, 'seasonal_yam')}

    def update_passage_categories(self, database, table, preview, index='accession', **kwargs):
        print("Updating passage_category field via passage field")
        upload = flu_upload(**args.__dict__)
        updated_sequences = []
        cursor = r.table(table).run()
        counter = 0
        total = r.table(table).count().run()
        print("Analyzing " +  str(total) + " sequence entries")
        for sequence in cursor:
            updated = upload.format_passage(sequence)
            if updated:
                updated_sequences.append(sequence)
            counter += 1
            if counter % 10000 == 0:
                print(str(counter) + " entries parsed")
        print("Found " +  str(len(updated_sequences)) + " sequences to update")
        if not preview:
            print("Updating " + str(len(updated_sequences)) + " sequence passage categories in " + database + "." + table)
            del kwargs['overwrite']
            self.upload_to_rethinkdb(database, table, updated_sequences, overwrite=True, index='accession')
        else:
            print("Preview of updates to be made, remove --preview to make updates to database")

    def update_groupings(self, viruses_table, sequences_table, database, update_keyword='tbd', preview=False, optimal_upload=50, **kwargs):
        '''
        Get viruses that have not had virus groupings determined, signaled by update_keyword
        Align HA sequences to outgroups to determine the closest grouping
        '''
        print("Updating grouping fields")
        print("Getting viruses from " + database + "." + viruses_table + " with grouping fields (vtype, subtype, or lineage) equal to ", update_keyword)
        viruses = list(r.db(database).table(viruses_table).filter((r.row["vtype"] == update_keyword) | (r.row["subtype"] == update_keyword) | (r.row["lineage"] == update_keyword)).run())
        ha_sequences = list(r.db(database).table(sequences_table).filter((r.row["locus"] == 'ha')).run())
        accession_to_sequence_doc = {doc['accession']:doc for doc in ha_sequences}
        # split updating of groups into groups of 50 viruses
        if len(viruses) > optimal_upload:
            list_viruses = [viruses[x:x+optimal_upload] for x in range(0, len(viruses), optimal_upload)]
        else:
            list_viruses = [viruses]
        print("Determining groupings for " + str(len(viruses)) + " viruses in " + str(len(list_viruses)) + " batches of " + str(optimal_upload) + " viruses at a time")
        for group_num, virus_group in enumerate(list_viruses):
            print("Group " + str(group_num+1)) + " out of " + str(len(list_viruses)) + " groupings"
            for virus in virus_group:
                for acc in virus['sequences']:
                    if acc in accession_to_sequence_doc:
                        # try determining grouping information from one HA sequence
                        result = self.align_flu(accession_to_sequence_doc[acc])
                        if result is not None:
                            virus['vtype'], virus['subtype'], virus['lineage'] = result
                        else:
                            virus['vtype'], virus['subtype'], virus['lineage'] = "undetermined", "undetermined", "undetermined"
                        break
            if not preview:
                print("Updating " + str(len(virus_group)) + " virus groupings in " + self.database + "." + self.viruses_table)
                self.upload_to_rethinkdb(self.database, self.viruses_table, virus_group, overwrite=True, optimal_upload=optimal_upload, index='strain')
            else:
                print("Preview of updates to be made, remove --preview to make updates to database")

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = flu_update(**args.__dict__)
    connVDB.update(**args.__dict__)
