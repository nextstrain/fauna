import os, re, time, datetime, csv, sys, json, errno
import pandas as pd
from upload import upload
from rethinkdb import r
from Bio import SeqIO
import argparse
import subprocess
from parse import parse
from upload import parser
import xlrd
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload
from titer_block import find_titer_block, find_serum_rows, find_virus_columns

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
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]

    if True in exten:
        ind = exten.index(True)
        convert_vidrl_xls_to_tsv(path, fstem, ind, assay_type)
    else:
        print("Unable to recognize file {}/{}".format(path,fstem))
        sys.exit()

def convert_vidrl_xls_to_tsv(path, fstem, ind, assay_type):
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])

    # Default patterns, VIDRL
    virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"
    serum_id_pattern = r"^[A-Z]\d{4,8}$"
    serum_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+)"
    serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"
    crick = False

    for worksheet in workbook.sheets():
        # autodetecting titer, strain, serum blocks
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

        serum_block = find_serum_rows(
            worksheet=worksheet,
            titer_coords=titer_coords,
            virus_names=virus_block["virus_names"],
            serum_id_pattern=serum_id_pattern,
            serum_passage_pattern=serum_passage_pattern,
            serum_abbrev_pattern=serum_abbrev_pattern,
            crick=crick,
        )

        # Print the most likely row and column indices for the titer block and the vote counts
        print("Titer block:")
        print(f"  Most likely (n={titer_block['col_start'][0][1]}) col_start: {titer_block['col_start'][0][0]}")
        print(f"  Most likely (n={titer_block['col_end'][0][1]}) col_end: {titer_block['col_end'][0][0]}")
        print(f"  Most likely (n={titer_block['row_start'][0][1]}) row_start: {titer_block['row_start'][0][0]}")
        print(f"  Most likely (n={titer_block['row_end'][0][1]}) row_end: {titer_block['row_end'][0][0]}")

        # Print virus and serum annotations row and column indices
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

        serum_mapping = serum_block["serum_mapping"]

        mat=worksheet

        with open('data/tmp/%s.tsv'%(fstem), 'w') as outfile:
            outfile.write("%s\n" % ("\t".join(ELIFE_COLUMNS)))

            print("assay_type: " + assay_type)
            # Zero-indexed positions
            row_start = titer_coords['row_start']
            row_end = titer_coords['row_end']
            col_start = titer_coords['col_start']
            col_end = titer_coords['col_end']

            virus_strain_col_index = virus_block['virus_col_idx']
            virus_passage_col_index = virus_block['virus_passage_col_idx']

            serum_id_row_index = serum_block['serum_id_row_idx']
            serum_passage_row_index = serum_block['serum_passage_row_idx']
            serum_strain_row_index = serum_block['serum_abbrev_row_idx']

            source = "vidrl_%s"%(fstem).strip()
            virus_passage_category = ''
            serum_passage_category = ''
            for i in range(row_start, (row_end+1)):
                virus_strain = str(mat.cell_value(i,virus_strain_col_index)).strip()
                virus_passage = str(mat.cell_value(i,virus_passage_col_index)).strip()
                for j in range(col_start, (col_end+1)):
                    serum_id = str(mat.cell_value(serum_id_row_index,j)).strip().replace(' ','')
                    serum_passage = str(mat.cell_value(serum_passage_row_index,j)).strip()
                    serum_abbr = str(mat.cell_value(serum_strain_row_index,j)).strip()
                    serum_abbr = serum_abbr.replace(' ','')
                    serum_strain = serum_mapping.get(serum_abbr, serum_abbr)
                    titer = str(mat.cell_value(i,j)).strip()
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
