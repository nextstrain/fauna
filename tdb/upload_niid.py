import os, re, time, datetime, csv, sys, json
from upload import upload
import rethinkdb as r
from Bio import SeqIO
import argparse
import subprocess
from parse import parse
from upload import parser
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload
import logging
logger = logging.getLogger()

def read_niid(path, fstem):
    '''
    Read all csv tables in path, create data frame with reference viruses as columns
    '''
    fname = path + fstem + ".csv"
    # import glob
    # flist = glob.glob(path + '/NIMR*csv') #BP
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]
    if True in exten:
        ind = exten.index(True)
        convert_xls_to_csv(path, fstem, ind)
        fname = "data/tmp/%s.csv"%(fstem)
        parse_niid_matrix_to_tsv(fname, path)
    else:
        logger.critical("Unable to recognize file extension of {}/{}".format(path,fstem))
        sys.exit()

def convert_xls_to_csv(path, fstem, ind):
    import xlrd
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for sheet in workbook.sheets():
        with open('data/tmp/%s.csv'%(fstem), 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(sheet.row_values(row) for row in range(10))
            temp = sheet.row_values(5)
            col = sheet.col_values(1)
            for i in range(4,12):
                temp[i] = col[i-3]
            writer.writerow(temp)
            writer.writerows(sheet.row_values(row) for row in range(11,sheet.nrows))
        return

def parse_niid_matrix_to_tsv(fname, original_path):
    from string import strip
    src_id = fname.split('/')[-1]
    with open(fname) as infile:
        csv_reader = csv.reader(infile)
        mat = list(csv_reader)
    with open('data/tmp/%s.tsv'%(src_id[:-4]), 'wb') as outfile:
        header = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]
        outfile.write("%s\n" % ("\t".join(header)))
        original_path = original_path.split('/')
        try:
            original_path.remove('')
        except:
            pass
        assay_type = determine_assay_type(original_path)
        for i in range(7,len(mat)):
            for j in range(4,12):
                virus_strain = mat[i][1]
                serum_strain = mat[5][j].split('\n')[0]
                serum_id = 'UnknownFerret'
                titer = mat[i][j]
                source = "niid_%s"%(src_id)
                virus_passage = mat[i][2]
                virus_passage_category = ''
                serum_passage = mat[5][j].split('\n')[1]
                serum_passage_category = ''
                line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                outfile.write(line)

def determine_assay_type(original_path):
    return original_path[-3]

def determine_subtype(original_path):
    original_path = original_path.split('/')
    try:
        original_path.remove('')
    except:
        pass
    subtype = original_path[-4]
    if subtype.lower() == "victoria":
        subtype = "vic"
    if subtype.upper() == "yamagata":
        subtype = "yam"
    return subtype

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    read_niid(args.path, args.fstem)
    subtype = determine_subtype(args.path)
    args.fstem = args.fstem.replace('(','\\(').replace(')','\\)')
    if args.preview:
        command = "python tdb/elife_upload.py -db niid_tdb --subtype " + subtype + " --path data/tmp/ --fstem " + args.fstem + " --preview"
        print command
        subprocess.call(command, shell=True)
    else:
        command = "python tdb/elife_upload.py -db niid_tdb --subtype " + subtype + " --path data/tmp/ --fstem " + args.fstem
        print command
        subprocess.call(command, shell=True)
