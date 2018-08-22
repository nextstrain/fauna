import os, re, time, datetime, csv, sys, json, errno
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
# logger = logging.getLogger()

parser.add_argument('--assay_type', default='hi')

def read_vidrl(path, fstem, assay_type):
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
        parse_vidrl_matrix_to_tsv(fname, path, assay_type)
    else:
        # logger.critical("Unable to recognize file extension of {}/{}".format(path,fstem))
        sys.exit()

def convert_xls_to_csv(path, fstem, ind):
    import xlrd
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for sheet in workbook.sheets():
        tmpfile = 'data/tmp/%s.csv'%(fstem)
        if not os.path.exists(os.path.dirname(tmpfile)):
            try:
                os.makedirs(os.path.dirname(tmpfile))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        with open(tmpfile, 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(sheet.row_values(row) for row in range(10))
            # Edit row containing serum strains
            # if '/' in sheet.row_values(10)[4].strip() and '/' in sheet.row_values(12)[2].strip():
            #     print sheet.row_values(12)[4]
            # else:
            #     raise ValueError
            # assume there all always 12 reference viruses in a table
            row_with_ref_sera = sheet.row_values(10)
            col_with_ref_viruses = sheet.col_values(2)
            ref_serum_indices = range(4, 16)
            ref_virus_indices = range(12, 24)
            for (sindex, vindex) in zip(ref_serum_indices, ref_virus_indices):
                row_with_ref_sera[sindex] = col_with_ref_viruses[vindex]
            writer.writerow(row_with_ref_sera)
            writer.writerows(sheet.row_values(row) for row in range(11,sheet.nrows))
        return

def parse_vidrl_matrix_to_tsv(fname, original_path, assay_type):
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
        print("assay_type: " + assay_type)
        if assay_type == "hi":
            start_row = 12
            start_col = 5 # Changed from 7 to 5 for VIC/YAM tables -BP
            end_col = 15
            virus_strain_col_index = 2
            virus_passage_col_index = 16
        elif assay_type == "fra":
            start_row = 12
            start_col = 4
            end_col = 13
            virus_strain_col_index = 2
            virus_passage_col_index = 14
            # some FRA tables have 10 sera, some have 11, some have 9
            check_cell_10th_sera = mat[start_col][13]
            check_cell_11th_sera = mat[start_col][14]
            if check_cell_10th_sera == '':
                virus_passage_col_index = 13
            elif check_cell_10th_sera != '' and check_cell_11th_sera == '':
                virus_passage_col_index = 14
            else:
                virus_passage_col_index = 15

        # some tables are do not begin where we think they do
        check_cell = mat[10][2]
        if check_cell == "Reference Antigens":
            for i in range(start_row, len(mat)):
                for j in range(start_col, end_col):
                    virus_strain = mat[i][virus_strain_col_index]
                    serum_strain = mat[10][j]
                    serum_id = mat[8][j]
                    titer = mat[i][j]
                    source = "vidrl_%s"%(src_id)
                    virus_passage = mat[i][virus_passage_col_index]
                    virus_passage_category = ''
                    serum_passage = mat[9][j]
                    serum_passage_category = ''
                    line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                    outfile.write(line)

# def determine_initial_indices(path, fstem):
#     import xlrd
#     exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]
#     if True in exten:
#         ind = exten.index(True)
#     else:
#         return
#     exts = ['.xls', '.xlsm', '.xlsx']
#     workbook = xlrd.open_workbook(path+fstem + exts[ind])
#     value = None
#     for sheet in workbook.sheets():
#         with open('data/tmp/%s.csv'%(fstem), 'wb') as f:
#             for row in range(len(sheet.col_values(0))):
#                 for col in range(len(sheet.row_values(0))):
#                     try:
#                         value = int(sheet.row_values(col)[row])
#                     except:
#                         pass
#                     if value:
#                         if is_power_of_2(value/10):
#                             return col, row
    # return 0, 0
#
# def is_power_of_2(n):
#     n = n/2
#     if n == 2:
#         return True
#     elif n > 2:
#         is_power_of_2(n)
#     else:
#         return False

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    else:
        if not args.path.endswith('/'):
            args.path = args.path + '/'
    if args.database is None:
        args.database = "vidrl_tdb"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    # x_shift, y_shift = determine_initial_indices(args.path, args.fstem)
    read_vidrl(args.path, args.fstem, args.assay_type)
    #TODO: This is where I will add conversion of vidrl files to eLife format!
    if args.subtype:
        if args.preview:
            command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + args.subtype + " --path data/tmp/ --fstem " + args.fstem + " --preview"
            print command
            subprocess.call(command, shell=True)
        else:
            command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + args.subtype + " --path data/tmp/ --fstem " + args.fstem
            print command
            subprocess.call(command, shell=True)
    else:
        print("Subtype needs to be specified with --subtype")
