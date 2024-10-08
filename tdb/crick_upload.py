import os, re, time, datetime, csv, sys, json
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

parser.add_argument('--assay_type', default='hi')

def build_location_mapping():
    l = { "Swit": "Switzerland",
          "Bris": "Brisbane",
          "Ire": "Ireland",
          "HK": "HongKong",
          "Maur": "Mauritius",
          "Nord-West": "NordrheinWestfalen",
          "Mich": "Michigan",
          "Bret": "Bretagne",
          "Catal": "Catalonia"}
    return l

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
            fname = "../fludata/Crick-London-WHO-CC/processed-data/csv/{}.csv".format(sheet)
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
        # comments sheets are just instructions on how to use the workbook template
        if sheet.name == 'comments':
            print(f"Skipping sheet {sheet.name!r}", file=sys.stderr)
            continue
        # Replace spaces with underscores in sheet name so that the call to
        # elife_upload does not error out due to the space in --fstem
        sheet.name = sheet.name.replace(' ', '_')
        with open('../fludata/Crick-London-WHO-CC/processed-data/csv/{}_{}.csv'.format(fstem, sheet.name), 'w') as f:
            writer = csv.writer(f)
            print(sheet.name)
            for row in range(sheet.nrows):
                new_row = []
                for cell in sheet.row_values(row):
                    try:
                        new_row.append(cell)
                    except:
                        import pdb; pdb.set_trace()
                writer.writerow(new_row)
        print("wrote new csv to ../fludata/Crick-London-WHO-CC/processed-data/csv/{}_{}.csv".format(fstem, sheet.name))
        sheets.append("{}_{}".format(fstem, sheet.name))
    return sheets

def parse_crick_matrix_to_tsv(fname, original_path, assay_type):
    src_id = fname.split('/')[-1]
    with open(fname) as infile:
        csv_reader = csv.reader(infile)
        mat = list(csv_reader)
    with open('../fludata/Crick-London-WHO-CC/processed-data/tsv/%s.tsv'%(src_id[:-4]), 'w') as outfile:
        header = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]
        outfile.write("%s\n" % ("\t".join(header)))
        original_path = original_path.split('/')
        try:
            original_path.remove('')
        except:
            pass
        if assay_type == "hi":
            start_row = 9
            start_col = 6
            col_span = 1
            virus_strain_col_index = 1
            virus_passage_col_index = 5
            serum_strain_row_index = 3
            serum_passage_row_index = 5
            serum_id_row_index = 6
        elif assay_type == "fra":
            start_row = 13
            start_col = 5
            col_span = 2
            virus_strain_col_index = 1
            virus_passage_col_index = 4
            serum_strain_row_index = 6
            serum_passage_row_index = 8
            serum_id_row_index = 9
        for i in range(start_row, len(mat)):
            for j in range(start_col, len(mat[0]), col_span):
                virus_strain = mat[i][virus_strain_col_index].strip()
                virus_strain = re.sub('\u0410', 'A', virus_strain) # Cyrillic A
                serum_strain = mat[serum_strain_row_index][j].rstrip("/")+"/"+mat[serum_strain_row_index+1][j].lstrip("/")
                m = build_location_mapping()
                for (k,v) in m.items():
                    if v not in serum_strain:
                        serum_strain = serum_strain.replace(k, v)
                serum_id = mat[serum_id_row_index][j]
                titer = mat[i][j]
                source = "crick_%s"%(src_id)
                virus_passage = mat[i][virus_passage_col_index]
                virus_passage_category = ''
                serum_passage = mat[serum_passage_row_index][j]
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
    if args.database is None:
        args.database = "crick_tdb"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    # x_shift, y_shift = determine_initial_indices(args.path, args.fstem)
    sheets = read_crick(args.path, args.fstem, args.assay_type)
    for sheet in sheets:
        if args.subtype:
            subtype = args.subtype
        else:
            subtype = determine_subtype(sheet)
        if args.preview:
            print("Subtype: {}".format(subtype))
            print("Sheet: {}".format(sheet))
            command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + subtype + " --path ../fludata/Crick-London-WHO-CC/processed-data/tsv/ --fstem " + sheet + " --preview"
            print(command)
            subprocess.call(command, shell=True)
        else:
            command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + subtype + " --path ../fludata/Crick-London-WHO-CC/processed-data/tsv/ --fstem " + sheet
            print(command)
            subprocess.call(command, shell=True)
