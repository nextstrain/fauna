import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import parser

class zibra_metadata_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def fix_name(self, name):
        return name

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = ['strain', 'amplicon_concentration', 'ct', 'date', 'division', 'lacen_id', 'location',
                    'microcephaly', 'onset_date', 'patient_age', 'patient_sex', 'rt_positive', 'sample_type']
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zibra_metadata_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
