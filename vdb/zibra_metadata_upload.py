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
    fasta_fields = ['strain', 'lacen_gal', 'rt_positive', 'ct', 'rt_date', 'amplicon_concentration', 'minion_barcodes',
                    'country', 'state', 'municipality', 'collection_date', 'onset_date', 'host_species', 'patient_age'
                    'patient_sex', 'pregnant', 'pregnancy_week', 'pregnancy_trimester', 'microcephaly', 'sample_type',
                    'symptoms']
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zibra_metadata_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
