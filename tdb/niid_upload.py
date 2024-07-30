import os, re, time, datetime, csv, sys, json
from upload import upload
from rethinkdb import r
from Bio import SeqIO
import argparse
import subprocess
import unicodedata
from parse import parse
import xlrd
from upload import parser
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload
from titer_block import find_titer_block, find_serum_rows, find_virus_columns

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
        convert_niid_xls_to_tsv(path, fstem, ind, subtype, assay_type)

def convert_niid_xls_to_tsv(path, fstem, ind, subtype, assay_type):
    # Set flutype
    suptype=subtype.lower()
    flutype = ""
    if subtype == "h3n2" or subtype == "h1n1pdm":
        flutype = "A"
    if subtype == "vic" or subtype == "yam":
        flutype = "B"

    # Set NIID patterns
    virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"
    serum_id_pattern = r".+(No\.|no\.).+"
    serum_passage_pattern = r".+(Egg|Cell).+"
    serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"
    crick = False

    # Open workbook
    wb_name = path + '/' + fstem + ind
    workbook = xlrd.open_workbook(filename=wb_name, encoding_override="cp1252")
    for worksheet_index, worksheet in enumerate(workbook.sheets(), start=1):
        print(f"Reading worksheet {worksheet_index} '{worksheet.name}' in file '{fstem}'")
        # autodetecting titer, virus, serum blocks
        titer_block = find_titer_block(worksheet)

        if len(titer_block["col_start"]) == 0:
            print("No titer block found.")
            break

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

        # If no virus names are found, might not be a valid worksheet, skip worksheet to avoid breaking find_serum_rows
        if virus_block["virus_names"] is None:
            print(f"Virus names not found. Check the virus pattern: '{virus_pattern}'")
            break

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
        print(f"Titer block: n = {titer_block['row_start'][0][1]}x{titer_block['col_start'][0][1]} = {titer_block['row_start'][0][1]*titer_block['col_start'][0][1]}")
        print(f"  Most likely (n={titer_block['col_start'][0][1]}) col_start: {titer_block['col_start'][0][0]}")
        print(f"  Most likely (n={titer_block['col_end'][0][1]}) col_end: {titer_block['col_end'][0][0]}")
        print(f"  Most likely (n={titer_block['row_start'][0][1]}) row_start: {titer_block['row_start'][0][0]}")
        print(f"  Most likely (n={titer_block['row_end'][0][1]}) row_end: {titer_block['row_end'][0][0]}")

        # For debugging purposes, print alternative indices (e.g. col_start, col_end, row_start, row_end)
        # print("Alternative indices:")
        # for i in range(1, len(titer_block['row_start'])):
        #     print(f"  Alternative (n={titer_block['row_start'][i][1]}) row_start: {titer_block['row_start'][i][0]}")

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

        serum_mapping = serum_block["serum_mapping"] # This is not used since serum_strain is being parsed in the loop below

        mat = worksheet

        with open('data/tmp/%s.tsv'%(fstem), 'w') as outfile:
            header = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]
            outfile.write("%s\n" % ("\t".join(header)))

            serum_id_row_index = serum_block['serum_id_row_idx']
            row_start = titer_coords['row_start']
            row_end = titer_coords['row_end']
            virus_id_col_index = virus_block['virus_col_idx']
            virus_passage_col_index=virus_block['virus_passage_col_idx']
            col_start = titer_coords['col_start']
            col_end = titer_coords['col_end']

            for i in range(row_start, row_end+1):
                for j in range(col_start, col_end+1):
                    virus_strain = str(mat.cell_value(i,virus_id_col_index)).strip()
                    serum_id = str(mat.cell_value(serum_id_row_index,j)).strip().replace(' ','')
                    serum_id = re.sub(r'[\r\n ]+', '', serum_id)
                    m = re.search(r'^(\S+)(egg|cell|siat|hck|nib121|ivr|\(bvr)', serum_id, re.IGNORECASE)
                    if m is None:
                        m = re.search(r'^(\S+)(no\.)', serum_id, re.IGNORECASE)
                    serum_strain = ""
                    if m:
                        serum_strain = m.group(1)
                        if not serum_strain.startswith(flutype + "/"):
                            serum_strain = flutype + "/" + serum_strain
                    # Normalize U+ff1c '＜' to U+003c '<'
                    titer = unicodedata.normalize('NFKC', str(mat.cell_value(i,j)).strip())
                    # Allow either "< 10" or "<10"
                    titer = re.sub(r'< ', '<', titer)
                    source = "niid_%s"%(fstem).strip()
                    virus_passage = str(mat.cell_value(i,virus_passage_col_index)).strip()
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
