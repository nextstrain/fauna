import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload

class Zika_vdb_upload(vdb_upload):

    def __init__(self, database, virus, source, fasta_fname, fasta_fields, locus):
        vdb_upload.__init__(self, database, virus, source, fasta_fname, fasta_fields, locus)

if __name__=="__main__":

    database = 'vdb'
    virus = 'Zika'
    source = 'virological'
    fasta_fname = 'zika_virological_02-22-2016.fasta'
    fasta_fields = {0:'accession', 1:'strain', 2:'date', 4:'country', 5:'division', 6:'location'}
    locus = 'Genome'


    upload_run = Zika_vdb_upload(database, virus, source, fasta_fname, fasta_fields, locus)
    upload_run.upload()
