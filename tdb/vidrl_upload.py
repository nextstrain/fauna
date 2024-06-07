import os, re, time, datetime, csv, sys, json, errno
import pandas as pd
from upload import upload
from rethinkdb import r
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

ELIFE_COLUMNS = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]

def parse_tsv_mapping_to_dict(tsv_file):
    map_dict = {}
    with open(tsv_file, 'r') as f:
        for line in f:
            (key, value) = line.split('\t')
            key = key.lower()
            map_dict[key] = value.rstrip('\n')
    return map_dict

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
        print("Unable to recognize file {}/{}".format(path,fstem))
        sys.exit()

def convert_xls_to_csv(path, fstem, ind):
    ref_sera_row=8    # Should be <= serum_strain_row_index
    titer_col_start=3 # Should be <= start_col
    import xlrd
    sera_mapping_file = 'source-data/vidrl_serum_mapping.tsv'
    sera_mapping = parse_tsv_mapping_to_dict(sera_mapping_file)
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
        with open(tmpfile, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(sheet.row_values(row) for row in range(10))
            # Edit row containing serum strains
            # if '/' in sheet.row_values(10)[4].strip() and '/' in sheet.row_values(12)[2].strip():
            #     print(sheet.row_values(12)[4])
            # else:
            #     raise ValueError
            # assume there all always 12 reference viruses in a table
            row_with_ref_sera = sheet.row_values(ref_sera_row)
            for i in range(titer_col_start,16):
                row_with_ref_sera[i] = row_with_ref_sera[i].strip('\n')
                try:
                    row_with_ref_sera[i] = sera_mapping[row_with_ref_sera[i].lower()]
                except KeyError:
                    print("Couldn't find {} in mapping lookup {}".format(row_with_ref_sera[i], sera_mapping_file))
                    with open ('data/BAD_VIDRL_KEYS.txt', 'a') as f:
                        f.write(row_with_ref_sera[i])
            writer.writerow(row_with_ref_sera)
            writer.writerows(sheet.row_values(row) for row in range(ref_sera_row+1,sheet.nrows))
        return

def parse_vidrl_matrix_to_tsv(fname, original_path, assay_type):
    src_id = fname.split('/')[-1]
    with open(fname) as infile:
        csv_reader = csv.reader(infile)
        mat = list(csv_reader)
    with open('data/tmp/%s.tsv'%(src_id[:-4]), 'w') as outfile:
        outfile.write("%s\n" % ("\t".join(ELIFE_COLUMNS)))
        original_path = original_path.split('/')
        try:
            original_path.remove('')
        except:
            pass
        print("assay_type: " + assay_type)
        if assay_type == "hi":
            # Zero-indexed positions
            start_row = 10
            start_col = 3
            end_col = 12+start_col # Changed from 7 to 5 for Vic/YAM tables -BP
            virus_strain_col_index = 1
            virus_passage_col_index = end_col
            serum_id_row_index = 6
            serum_passage_row_index = 7
            serum_strain_row_index = 8
        elif assay_type == "fra":
            start_row = 12
            start_col = 4
            end_col = 13
            virus_strain_col_index = 2
            virus_passage_col_index = 14
            serum_id_row_index = 8
            serum_passage_row_index = 9
            serum_strain_row_index = 10
            # # some FRA tables have 10 sera, some have 11, some have 9
            # check_cell_10th_sera = mat[start_col][13]
            # check_cell_11th_sera = mat[start_col][14]
            # if check_cell_10th_sera == '':
            #     virus_passage_col_index = 13
            # elif check_cell_10th_sera != '' and check_cell_11th_sera == '':
            #     virus_passage_col_index = 14
            # else:
            #     virus_passage_col_index = 15

        # some tables are do not begin where we think they do
        # add possible starting locations
        possible_starts = [mat[serum_strain_row_index][virus_strain_col_index], mat[10][2], mat[13][0]]
        for check_cell in possible_starts:
            if check_cell == "Reference Antigens":
                for i in range(start_row, len(mat)):
                    for j in range(start_col, end_col):
                        virus_strain = mat[i][virus_strain_col_index].strip()
                        serum_strain = mat[serum_strain_row_index][j].strip()
                        serum_id = mat[serum_id_row_index][j].strip().replace(' ','')
                        titer = mat[i][j].strip()
                        source = "vidrl_%s"%(src_id).strip()
                        virus_passage = mat[i][virus_passage_col_index].strip()
                        virus_passage_category = ''
                        serum_passage = mat[serum_passage_row_index][j].strip()
                        serum_passage_category = ''
                        line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                        outfile.write(line)


def read_flat_vidrl(path, fstem, assay_type):
    """
    Read the flat CSV file with *fstem* in the provided *path* and convert
    to the expected TSV file at `data/tmp/<fstem>.tsv` for tdb/elife_upload.
    """
    column_map = parse_tsv_mapping_to_dict("source-data/vidrl_flat_file_column_map.tsv")
    filepath = path + fstem + ".csv"

    titer_measurements = pd.read_csv(filepath, usecols=column_map.keys()) \
                           .rename(columns=column_map)

    titer_measurements["assay_type"] = assay_type
    titer_measurements["virus_passage_category"] = ""
    titer_measurements["serum_passage_category"] = ""
    titer_measurements["source"] = "vidrl_{}.csv".format(fstem)

    titer_measurements[ELIFE_COLUMNS].to_csv("data/tmp/{}.tsv".format(fstem), sep="\t", index=False)


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

    if args.ftype == "flat":
        read_flat_vidrl(args.path, args.fstem, args.assay_type)
    else:
        read_vidrl(args.path, args.fstem, args.assay_type)

    if args.subtype:
        if args.preview:
            command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + args.subtype + " --path data/tmp/ --fstem " + args.fstem + " --preview"
            print(command)
            subprocess.call(command, shell=True)
        else:
            command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + args.subtype + " --path data/tmp/ --fstem " + args.fstem
            print(command)
            subprocess.call(command, shell=True)
    else:
        print("Subtype needs to be specified with --subtype")
