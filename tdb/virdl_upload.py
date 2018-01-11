import os, re, time, datetime, csv, sys, json
from upload import upload
import rethinkdb as r
from Bio import SeqIO
import argparse
from parse import parse
from upload import parser
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload

class nimr_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.assay_type = kwargs['assay_type']
        self.assay_date = None

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = virdl_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
