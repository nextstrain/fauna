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

def parse_sera_mapping_to_dict():
    sera_mapping_file = 'source-data/vidrl_serum_mapping.tsv'
    map_dict = {}
    with open(sera_mapping_file, 'r') as f:
        for line in f:
            (key, value) = line.split('\t')
            map_dict[key] = value
    return map_dict

def read_vidrl(path, fstem):
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
        parse_vidrl_matrix_to_tsv(fname, path)
    else:
        logger.critical("Unable to recognize file extension of {}/{}".format(path,fstem))
        sys.exit()

def convert_xls_to_csv(path, fstem, ind):
    import xlrd

    sera_mapping = parse_sera_mapping_to_dict()

    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for sheet in workbook.sheets():
        with open('data/tmp/%s.csv'%(fstem), 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(sheet.row_values(row) for row in range(10))
            # Edit row containing serum strains
            print sheet.row_values(10)[4]
            print type(sheet.row_values(10)[4])
            print sheet.row_values(12)[2]
            print type(sheet.row_values(12)[2])
            if '/' in sheet.row_values(10)[4].strip() and '/' in sheet.row_values(12)[2].strip():
                print sheet.row_values(12)[4]
            else:
                raise ValueError
            ### This is the bit that needs to be replaced
            temp = sheet.row_values(10)
            for i in range(4,16):
                temp[i] = temp[i].strip('\n')
                try:
                    temp[i] = sera_mapping[temp[i]]
                except KeyError:
                    pass
            writer.writerow(temp)
            ###
            writer.writerows(sheet.row_values(row) for row in range(11,sheet.nrows))
        return

def parse_vidrl_matrix_to_tsv(fname, original_path):
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

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    # x_shift, y_shift = determine_initial_indices(args.path, args.fstem)
    read_vidrl(args.path, args.fstem)
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
