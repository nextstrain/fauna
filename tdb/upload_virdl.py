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

def read_virdl(path, fstem, **kwargs):
    '''
    Read all csv tables in path, create data frame with reference viruses as columns
    '''
    fname = path + fstem + ".csv"
    # import glob
    # flist = glob.glob(path + '/NIMR*csv') #BP
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', 'xlsx'] ]
    if True in exten:
        ind = exten.index(True)
        convert_xls_to_csv(path, fstem, ind)
        fname = "data/tmp/%s.csv"%(fstem)
        parse_virdl_matrix_to_tsv(fname, path)

def convert_xls_to_csv(path, fstem, ind):
    import xlrd
    exts = ['.xls', '.xlsm', 'xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for sheet in workbook.sheets():
        with open('data/tmp/%s.csv'%(fstem), 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(sheet.row_values(row) for row in range(sheet.nrows))

def parse_virdl_matrix_to_tsv(fname, original_path):
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
                source = "VIRDL_%s"%(src_id)
                virus_passage = mat[i][16]
                virus_passage_category = ''
                serum_passage = mat[9][j]
                serum_passage_category = ''
                line = "%s\n" % ("\t".join([ virus_strain, serum_strain, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                outfile.write(line)

if __name__=="__main__":
    args = parser.parse_args()
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    read_virdl(args.path, args.fstem)
    print "Done with that step!"
    #TODO: This is where I will add conversion of VIRDL files to eLife format!
    sys.exit()
    connTDB = elife_upload(**args.__dict__)
    connTDB.upload(**args.__dict__)
