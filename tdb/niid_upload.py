import os, re, time, datetime, csv, sys, json
from upload import upload
from rethinkdb import r
from Bio import SeqIO
import argparse
import subprocess
import unicodedata
from parse import parse
from upload import parser
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload
import logging
# logger = logging.getLogger()

parser.add_argument('--assay_type', default='hi')

def read_niid(path, fstem, subtype, assay_type):
    '''
    Convert xls tables to csv tables, then parse to flat tsv
    '''
    possible_files = [ path + '/' + fstem + ext for ext in ['.xls', '.xlsm', '.xlsx']]
    real_file = ''
    for possible_file in possible_files:
        if os.path.isfile(possible_file):
            real_file = possible_file
    if real_file != '':
        print("real_file: " + real_file)
        ind = '.{}'.format(real_file.split('.')[-1])
        convert_xls_to_csv(path, fstem, ind)
        fname = "data/tmp/{}.csv".format(fstem)
        parse_niid_matrix_to_tsv(fname, path, subtype, assay_type)

def convert_xls_to_csv(path, fstem, ind):
    import xlrd
    wb_name = path + '/' + fstem + ind
    workbook = xlrd.open_workbook(filename=wb_name, encoding_override="cp1252")
    for sheet in workbook.sheets():
        with open('data/tmp/%s.csv'%(fstem), 'w') as f:
            writer = csv.writer(f)
            rows = []
            for i in range(sheet.nrows):
                row = []
                for j in range(sheet.ncols):
                    val = sheet.cell_value(i, j)
                    row.append(val)
                rows.append(row)
            writer.writerows(rows)
        return

def parse_niid_matrix_to_tsv(fname, original_path, subtype, assay_type):
    suptype=subtype.lower()
    flutype = ""
    if subtype == "h3n2" or subtype == "h1n1pdm":
        flutype = "A"
    if subtype == "vic" or subtype == "yam":
        flutype = "B"
    src_id = fname.split('/')[-1]
    with open(fname) as infile:
        csv_reader = csv.reader(infile)
        mat = list(csv_reader)
    with open('data/tmp/%s.tsv'%(src_id[:-4]), 'w') as outfile:
        header = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]
        outfile.write("%s\n" % ("\t".join(header)))
        original_path = original_path.split('/')
        try:
            original_path.remove('')
        except:
            pass
        if subtype == "h3n2":
            start_row = 7
            start_col = 4
            serum_id_row_index = 5
        elif subtype == "h1n1pdm":
            start_row = 6
            start_col = 4
            serum_id_row_index = 4
        elif subtype == "vic":
            start_row = 5
            start_col = 4
            serum_id_row_index = 4
        elif subtype == "yam":
            start_row = 5
            start_col = 4
            serum_id_row_index = 3
        for i in range(start_row, len(mat)):
            for j in range(start_col, len(mat[0])):
                virus_strain = mat[i][1]
                serum_id = mat[serum_id_row_index][j]
                serum_id = re.sub(r'[\r\n ]+', '', serum_id)
                m = re.search(r'^(\S+)(egg|cell|siat|hck|nib121|ivr|\(bvr)', serum_id, re.IGNORECASE)
                if m is None:
                    m = re.search(r'^(\S+)(no\.)', serum_id, re.IGNORECASE)
                serum_strain = ""
                # import pdb; pdb.set_trace()
                if m:
                    serum_strain = m.group(1)
                    if not serum_strain.startswith(flutype + "/"):
                        serum_strain = flutype + "/" + serum_strain
                # Normalize U+ff1c '＜' to U+003c '<'
                titer = unicodedata.normalize('NFKC', mat[i][j])
                source = "niid_%s"%(src_id)
                virus_passage = mat[i][2]
                virus_passage_category = ''
                serum_passage = "unknown"
                m = re.search(r'(egg)', serum_id, re.IGNORECASE)
                if m:
                    serum_passage = m.group(1)
                m = re.search(r'(cell|siat|hck)', serum_id, re.IGNORECASE)
                if m:
                    serum_passage = m.group(1)
                serum_passage_category = ''
                line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                outfile.write(line)

def determine_subtype(original_path):
    original_path = original_path.lower().split('/')
    if 'h3n2' in original_path:
        subtype = 'h3n2'
    elif 'h1n1pdm' in original_path:
        subtype = 'h1n1pdm'
    elif 'victoria' in original_path:
        subtype = 'vic'
    elif 'yamagata' in original_path:
        subtype = 'yam'
    else:
        subtype = "UnknownSubtype"
    return subtype

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if args.database is None:
        args.database = "niid_tdb"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    subtype = determine_subtype(args.path)
    read_niid(args.path, args.fstem, subtype, args.assay_type)
    args.fstem = args.fstem.replace('(','\\(').replace(')','\\)')
    if args.preview:
        command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + subtype + " --path data/tmp/ --fstem " + args.fstem + " --preview"
        print(command)
        subprocess.call(command, shell=True)
    else:
        command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + subtype + " --path data/tmp/ --fstem " + args.fstem
        print(command)
        subprocess.call(command, shell=True)
