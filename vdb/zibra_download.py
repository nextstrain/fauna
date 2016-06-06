import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from download import download
from download import parser

class zibra_download(download):
    def __init__(self, **kwargs):
        download.__init__(self, **kwargs)
        self.virus_specific_fasta_fields = []

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = ['strain', 'amplicon_concentration', 'ct', 'date', 'division', 'location', 'microcephaly',
                    'onset_date', 'patient_age', 'patient_sex', 'rt_positive', 'sample_type']
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zibra_download(**args.__dict__)
    connVDB.download(**args.__dict__)
