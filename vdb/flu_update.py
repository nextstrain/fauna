import os
import numpy as np
import rethinkdb as r
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

    def update_groupings(self, viruses_table, sequences_table, database, preview, **kwargs):
        print("Updating grouping fields")
        viruses = list(r.db(database).table(viruses_table).filter((r.row["vtype"] == 'tbd') | (r.row["subtype"] == 'tbd') | (r.row["lineage"] == 'tbd')).run())
        ha_sequences = list(r.db(database).table(sequences_table).filter((r.row["locus"] == 'ha')).run())
        accession_to_sequence_doc = {doc['accession']:doc for doc in ha_sequences}
        # split updating of groups into groups of 50 viruses
        optimal_upload = 50
        if len(viruses) > optimal_upload:
            list_viruses = [viruses[x:x+optimal_upload] for x in range(0, len(viruses), optimal_upload)]
        else:
            list_viruses = [viruses]
        print("Determining groupings for " + str(len(viruses)) + " viruses in " + str(len(list_viruses)) + " batches")
        for group in list_viruses:
            print("Determining groupings for " + str(len(group)) + " viruses")
            for virus in group:
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
                print("Updating " + str(len(group)) + " virus groupings in " + self.database + "." + self.viruses_table)
                self.upload_to_rethinkdb(self.database, self.viruses_table, group, 'update')
            else:
                print("Preview of updates to be made, remove --preview to make updates to database")

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = flu_update(**args.__dict__)
    connVDB.update(**args.__dict__)
