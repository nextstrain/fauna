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
# import logging
# print 'yay'
# logger = logging.getLogger()
# print 'more yay'

parser.add_argument('--assay_type', default='hi')

def read_crick(path, fstem, assay_type):
    '''
    Read all csv tables in path, create data frame with reference viruses as columns
    '''
    fname = path + fstem # + ".csv"
    # import glob
    # flist = glob.glob(path + '/NIMR*csv') #BP
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]
    if True in exten:
        ind = exten.index(True)
        sheets = convert_xls_to_csv(path, fstem, ind)
        for sheet in sheets:
            fname = "data/tmp/{}.csv".format(sheet)
            parse_crick_matrix_to_tsv(fname, path, assay_type)
    else:
        # logger.critical("Unable to recognize file extension of {}/{}".format(path,fstem))
        print("EXITING")
        sys.exit()
    return sheets

def convert_xls_to_csv(path, fstem, ind):
    import xlrd
    sheets = []
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for sheet in workbook.sheets():
        with open('data/tmp/{}_{}.csv'.format(fstem, sheet.name), 'wb') as f:
            try:
                writer = csv.writer(f)
                print(sheet.name)
                writer.writerows( sheet.row_values(row) for row in range(sheet.nrows))
                # writer.writerow(temp)
                # writer.writerows(sheet.row_values(row) for row in range(11,sheet.nrows))
            except:
                print("couldn't write data/tmp/{}_{}.csv".format(fstem,sheet.name))
        print("wrote new csv to data/tmp/{}_{}.csv".format(fstem, sheet.name))
        sheets.append("{}_{}".format(fstem, sheet.name))
        # sys.exit()
    return sheets

def parse_crick_matrix_to_tsv(fname, original_path, assay_type):
    assay_type = assay_type
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
        for i in range(14,len(mat)):
            for j in range(7,len(mat[0])):
                virus_strain = mat[i][1]
                serum_strain = mat[6][j]
                serum_id = mat[9][j]
                titer = mat[i][j]
                source = "crick_%s"%(src_id)
                virus_passage = mat[i][6]
                virus_passage_category = ''
                serum_passage = mat[8][j]
                serum_passage_category = ''
                line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                outfile.write(line)

def determine_subtype(fname):
    if fname.lower().startswith('h3n2'):
        subtype = 'h3n2'
    elif fname.lower().startswith ('h1n1pdm'):
        subtype = 'h1n1pdm'
    elif fname.lower().startswith('bvic'):
        subtype = 'vic'
    elif fname.lower().startswith('byam'):
        subtype = 'yam'
    else:
        subtype = 'unknown'
    return subtype


if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    # x_shift, y_shift = determine_initial_indices(args.path, args.fstem)
    sheets = read_crick(args.path, args.fstem, args.assay_type)
    #TODO: This is where I will add conversion of crick files to eLife format!
    for sheet in sheets:
        subtype = determine_subtype(sheet)
        if args.preview:
            print("Subtype: {}".format(subtype))
            print("Sheet: {}".format(sheet))
            command = "python tdb/elife_upload.py -db crick_tdb --subtype " + subtype + " --path data/tmp/ --fstem " + sheet + " --preview"
            print command
            subprocess.call(command, shell=True)
        else:
            command = "python tdb/elife_upload.py -db crick_tdb --subtype " + subtype + " --path data/tmp/ --fstem " + sheet
            print command
            subprocess.call(command, shell=True)
