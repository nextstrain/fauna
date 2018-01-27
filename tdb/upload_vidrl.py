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

def read_vidrl(path, fstem, x_shift, y_shift, **kwargs):
    '''
    Read all csv tables in path, create data frame with reference viruses as columns
    '''
    fname = path + fstem + ".csv"
    # import glob
    # flist = glob.glob(path + '/NIMR*csv') #BP
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]
    if True in exten:
        ind = exten.index(True)
        convert_xls_to_csv(path, fstem, ind, x_shift, y_shift)
        fname = "data/tmp/%s.csv"%(fstem)
        parse_vidrl_matrix_to_tsv(fname, path)
    else:
        logger.critical("Unable to recognize file extension of {}/{}".format(path,fstem))
        sys.exit()

def convert_xls_to_csv(path, fstem, ind, x_shift, y_shift):
    import xlrd
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for sheet in workbook.sheets():
        with open('data/tmp/%s.csv'%(fstem), 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(sheet.row_values(row) for row in range(10))
            # Edit row containing serum strains
            temp = sheet.row_values(y_shift-2)
            col = sheet.col_values(x_shift-2)
            # print col
            # print x_shift
            # print y_shift
            # print '#####'
            # print sheet.col_values(x_shift-1)
            # print '#####'
            # print sheet.col_values(x_shift-3)
            # sys.exit()
            for i in range(4,16):
                temp[i] = col[11+i-4]
            writer.writerow(temp)
            writer.writerows(sheet.row_values(row) for row in range(11,sheet.nrows))

def parse_vidrl_matrix_to_tsv(fname, original_path, x_shift, y_shift):
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
        assay_type = original_path[-1]
        for i in range(12,len(mat)):
            for j in range(4,15):
                virus_strain = mat[i][2]
                serum_strain = mat[10][j]
                serum_id = mat[8][j]
                titer = mat[i][j]
                source = "vidrl_%s"%(src_id)
                virus_passage = mat[i][16]
                virus_passage_category = ''
                serum_passage = mat[9][j]
                serum_passage_category = ''
                line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                outfile.write(line)

def determine_subtype(original_path):
    original_path = original_path.split('/')
    try:
        original_path.remove('')
    except:
        pass
    subtype = original_path[-2]
    return subtype

def determine_initial_indices(path, fstem):
    import xlrd
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]
    if True in exten:
        ind = exten.index(True)
    else:
        return
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    value = None
    for sheet in workbook.sheets():
        with open('data/tmp/%s.csv'%(fstem), 'wb') as f:
            for row in range(len(sheet.col_values(0))):
                for col in range(len(sheet.row_values(0))):
                    try:
                        value = int(sheet.row_values(col)[row])
                    except:
                        pass
                    if value:
                        if is_power_of_2(value/10):
                            return col, row
    return 0, 0

def is_power_of_2(n):
    n = n/2
    if n == 2:
        return True
    elif n > 2:
        is_power_of_2(n)
    else:
        return False

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    x_shift, y_shift = determine_initial_indices(args.path, args.fstem)
    read_vidrl(args.path, args.fstem, x_shift, y_shift)
    ####
    subtype = determine_subtype(args.path)
    #TODO: This is where I will add conversion of vidrl files to eLife format!
    if args.preview:
        command = "python tdb/elife_upload.py -db vidrl_tdb --subtype " + subtype + " --path data/tmp/ --fstem " + args.fstem + " --preview"
        print command
        subprocess.call(command, shell=True)
    else:
        command = "python tdb/elife_upload.py -db vidrl_tdb --subtype " + subtype + " --path data/tmp/ --fstem " + args.fstem
        print command
        subprocess.call(command, shell=True)
