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
from titer_block import find_titer_block, find_serum_rows, find_virus_columns

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
        sheets=convert_crick_xls_to_tsv(path, fstem, ind, assay_type)

        #sheets = convert_xls_to_csv(path, fstem, ind)
        #for sheet in sheets:
        #    fname = "../fludata/Crick-London-WHO-CC/processed-data/csv/{}.csv".format(sheet)
        #    parse_crick_matrix_to_tsv(fname, path, assay_type)
    else:
        # logger.critical("Unable to recognize file extension of {}/{}".format(path,fstem))
        print("EXITING")
        sys.exit()
    return sheets

def convert_crick_xls_to_tsv(path, fstem, ind, assay_type):
    # Set Crick patterns
    virus_pattern = r"[A-Z]/[\w\s-]+"
    virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"
    serum_id_pattern = r"F\d+/\d+"
    serum_passage_pattern = r"(MDCK|SIAT|Egg)"
    serum_abbrev_pattern = r"[A-Z]/[\w\s-]+"
    crick = True

    import xlrd
    sheets = []
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])
    for worksheet in workbook.sheets():
        worksheet.name = worksheet.name.replace(' ', '_')
        # autodetecting titer, strain, serum blocks
        titer_block = find_titer_block(worksheet)

        if len(titer_block["col_start"]) == 0 and assay_type == "hi":
            print("No titer block found.")
            break
        else:
            sheets.append("{}_{}".format(fstem, worksheet.name))

        titer_coords = {
            'col_start': titer_block["col_start"][0][0],
            'col_end': titer_block["col_end"][0][0],
            'row_start': titer_block["row_start"][0][0],
            'row_end': titer_block["row_end"][0][0]
        }

        virus_block = find_virus_columns(
            worksheet=worksheet,
            titer_coords=titer_coords,
            virus_pattern=virus_pattern,
            virus_passage_pattern=virus_passage_pattern,
        )

        serum_block = find_serum_rows(
            worksheet=worksheet,
            titer_coords=titer_coords,
            virus_names=virus_block["virus_names"],
            serum_id_pattern=serum_id_pattern,
            serum_passage_pattern=serum_passage_pattern,
            serum_abbrev_pattern=serum_abbrev_pattern,
            crick=crick,
        )

        # Print the most likely row and column indices for the titer block
        print("Titer block:")
        print(f"  Most likely (n={titer_block['col_start'][0][1]}) col_start: {titer_block['col_start'][0][0]}")
        print(f"  Most likely (n={titer_block['col_end'][0][1]}) col_end: {titer_block['col_end'][0][0]}")
        print(f"  Most likely (n={titer_block['row_start'][0][1]}) row_start: {titer_block['row_start'][0][0]}")
        print(f"  Most likely (n={titer_block['row_end'][0][1]}) row_end: {titer_block['row_end'][0][0]}")

        # Print virus and serum annotations row and column indices
        # Print Virus and Serum annotations row and column indices
        print("Virus (antigen) block: left and right of the titer block")
        print(f"  virus column index: {virus_block['virus_col_idx']}")
        print(f"  virus passage column index: {virus_block['virus_passage_col_idx']}")
        print(f"  virus names: {virus_block['virus_names']}")

        print("Serum (antisera) block: above the titer block")
        print(f"  serum ID row index: {serum_block['serum_id_row_idx']}")
        print(f"  serum passage row index: {serum_block['serum_passage_row_idx']}")
        print(f"  serum abbreviated name row index: {serum_block['serum_abbrev_row_idx']}")

        # Match abbreviated names across the top to the full names along the left side and auto convert to full names
        if serum_block["serum_abbrev_row_idx"] is not None:
            print("serum_mapping = {")
            for abbrev, full in serum_block["serum_mapping"].items():
                print(f"    '{abbrev}': '{full}',")
            print("}")

        serum_mapping = serum_block["serum_mapping"] # unused

        mat=worksheet

        if assay_type == "hi":
            start_row=titer_coords['row_start']
            start_col=titer_coords['col_start']
            end_row=titer_coords['row_end']
            end_col=titer_coords['col_end']
            virus_strain_col_index = virus_block['virus_col_idx']
            virus_passage_col_index = virus_block['virus_passage_col_idx']
            serum_strain_row_index = serum_block['serum_abbrev_row_idx']
            serum_passage_row_index = serum_block['serum_passage_row_idx']
            serum_id_row_index = serum_block['serum_id_row_idx']
            col_span = 1
        elif assay_type == "fra":
            start_row = 13
            start_col = 5
            end_row = len(mat)
            end_col = len(mat[0])
            col_span = 2
            virus_strain_col_index = 1
            virus_passage_col_index = 4
            serum_strain_row_index = 6
            serum_passage_row_index = 8
            serum_id_row_index = 9

        with open('../fludata/Crick-London-WHO-CC/processed-data/tsv/%s.tsv'%(fstem), 'w') as outfile:
            header = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]
            outfile.write("%s\n" % ("\t".join(header)))

            for i in range(start_row, end_row+1):
                for j in range(start_col, end_col+1, col_span):
                    virus_strain = str(mat.cell_value(i,virus_strain_col_index)).strip()
                virus_strain = re.sub('\u0410', 'A', virus_strain) # Cyrillic A
                    serum_strain = str(mat.cell_value(serum_strain_row_index,j)).strip().rstrip("/")+"/"+str(mat.cell_value(serum_strain_row_index+1,j)).strip().lstrip("/")
                    m = build_location_mapping()
                    for (k,v) in m.items():
                        if v not in serum_strain:
                            serum_strain = serum_strain.replace(k, v)
                    serum_id = str(mat.cell_value(serum_id_row_index,j)).strip()
                    titer = str(mat.cell_value(i,j)).strip()
                    source = "crick_%s"%(fstem).strip()
                    virus_passage = str(mat.cell_value(i,virus_passage_col_index)).strip()
                    virus_passage_category = ''
                    serum_passage = str(mat.cell_value(serum_passage_row_index,j)).strip()
                    serum_passage_category = ''
                    line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                    outfile.write(line)

    return(sheets)

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
