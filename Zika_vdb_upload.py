import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from vdb_upload import vdb_upload
from vdb_upload import parser

class Zika_vdb_upload(vdb_upload):

    def __init__(self,  fasta_fields, fasta_fname, database, virus, source, locus=None, vsubtype=None, authors=None, path=None, auth_key=None):
        '''
        :param fasta_fields: Dictionary defining position in fasta field to be included in database
        '''
        vdb_upload.__init__(self, fasta_fields, fasta_fname, database, virus, source, locus, vsubtype, authors, path, auth_key)

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 2:'strain', 4:'date', 6:'country'}
    #  0        1          2                            3  4       5     6
    # >KU647676|Zika_virus|MRS_OPY_Martinique_PaRi_2015|NA|2015_12|Human|Martinique
    run = Zika_vdb_upload(fasta_fields, fasta_fname=args.fname, database=args.database, virus=args.virus, source=args.source, locus=args.locus, vsubtype=args.subtype, authors=args.authors, path=args.path)
    run.upload()