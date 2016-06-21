import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from download import download
from download import get_parser

class zibra_download(download):
    def __init__(self, **kwargs):
        download.__init__(self, **kwargs)
        self.virus_specific_fasta_fields = []

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    fasta_fields = ['strain', 'lacen_gal', 'rt_positive', 'ct', 'rt_date', 'amplicon_concentration', 'minion_barcodes',
                    'country', 'state', 'municipality', 'collection_date', 'onset_date', 'host_species', 'patient_age'
                    'patient_sex', 'pregnant', 'pregnancy_week', 'pregnancy_trimester', 'microcephaly', 'sample_type',
                    'symptoms']
    setattr(args, 'fasta_fields', fasta_fields)
    connVDB = zibra_download(**args.__dict__)
    connVDB.download(**args.__dict__)
