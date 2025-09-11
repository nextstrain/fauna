#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import argparse
import xlrd
from collections import defaultdict
import re
sys.path.append('')  # need to import from base
from titer_block import find_titer_block, find_serum_rows, find_virus_columns

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Flatten a block of titers in an Excel worksheet.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--file",
        default="data/today/12-6-24\ Egg\ REF\ RK.xlsx",
        required=False,
        help="Path to the Excel file",
    )

    return parser.parse_args()

def main():
    args = parse_args()

    # Default patterns, VIDRL
    virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"
    serum_id_pattern = r"^[A-Z]\d{4,8}"
    serum_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+|CELL|cell|Cell)"
    serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"
    crick = False

    # Load the Excel file
    workbook = xlrd.open_workbook(args.file)
    filebasename = os.path.splitext(os.path.basename(args.file))[0]

    for worksheet_index, worksheet in enumerate(workbook.sheets(), start=1):
        print(f"Reading worksheet {worksheet_index} '{worksheet.name}' in file '{args.file}'")

        # Find the block of titers in the worksheet
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
            log_human_sera=False
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
            # print("Serum mapping:")
            # for abbrev, full in serum_block["serum_mapping"].items():
            #     print(f"  {abbrev} -> {full}")

            print("serum_mapping = {")
            for abbrev, full in serum_block["serum_mapping"].items():
                print(f"    '{abbrev}': '{full}',")
            print("}")


        # Check if all the necessary indices were found
        if virus_block["virus_col_idx"] is None:
            print(f"Virus column index not found. Check the virus pattern: '{virus_pattern}'")

        if virus_block["virus_passage_col_idx"] is None:
            print(f"Virus passage column index not found. Check the virus passage pattern: '{virus_passage_pattern}'")

        if serum_block["serum_id_row_idx"] is None:
            print(f"Serum ID row index not found. Check the serum ID pattern: '{serum_id_pattern}'")

        if serum_block["serum_passage_row_idx"] is None:
            print(f"Serum passage row index not found. Check the serum passage pattern: '{serum_passage_pattern}'")

        if serum_block["serum_abbrev_row_idx"] is None:
            print(f"Serum abbreviated name row index not found. Check the serum abbreviated name pattern: '{serum_abbrev_pattern}'")

        # flatten to tsv
        serum_mapping = serum_block["serum_mapping"]

        with open(f"{filebasename.replace(' ', '_')}_{worksheet.name.replace(' ', '_')}.tsv", 'w') as outfile:
          outfile.write("%s\n" % ("\t".join(["ag_strain_name","sr_strain_name","titer_value"])))
          # print
          for i in range(titer_coords['row_start'], titer_coords['row_end']):
              ag_strain_name=str(worksheet.cell_value(i,virus_block['virus_col_idx'])).strip()
              for j in range(titer_coords['col_start'], titer_coords['col_end']):
                  sr_strain_abbr=str(worksheet.cell_value(serum_block['serum_abbrev_row_idx'],j)).strip()
                  sr_strain_name=serum_mapping.get(sr_strain_abbr, sr_strain_abbr)
                  titer = str(worksheet.cell_value(i,j)).strip()
                  outfile.write(f"{ag_strain_name}\t{sr_strain_name}\t{titer}\n")

if __name__ == "__main__":
    main()
