#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys
import argparse
import xlrd
from datetime import datetime
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
    parser.add_argument(
        "--test-protocol",
        default="hi_protocol",
        required=False,
        help="Pass in assay protocol",
    )
    parser.add_argument(
        "--lot-column",
        default=3,
        required=False,
        help="Numeric column which has the lot id, below the titer block",
    )
    parser.add_argument(
        "--subtype-column",
        default=22,
        required=False,
        help="Subtype column, since sometimes to the left or right of the titer block or even intersperced",
        # TODO: autodetect this later
    )
    parser.add_argument(
        "--cdcid-column",
        default=19,
        required=False,
        help="CDCID column, since sometimes to the left or right of the titer block",
        # TODO: Autodetect this later
    )
    parser.add_argument(
        "--test-date",
        default="YYYY-MM-DD",
        required=False,
        help="Default test date, tends to be reported in very inconsistent places in the excel file",
    )
    parser.add_argument(
        "--col-end-adjustment",
        default=-3,
        required=False,
        help="Since a pair of numeric columns occur after the titer block, adjustment of the end of the titer block",
    )
    return parser.parse_args()

def is_titer_value(value):
    """
    Check if the value is numeric or a string representing a numeric value with '<'.

    Basically checking for titer values e.g. "80", "< 10", "1400", ">2560", etc.
    """
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, str):
        value = value.strip()
        if value.startswith("<") or value.startswith('＜'):
            try:
                float(value[1:].strip())
                return True
            except ValueError:
                return False
        if value.startswith(">"):
            try:
                float(value[1:].strip())
                return True
            except ValueError:
                return False
    return False

def find_sr_lookup(worksheet, titer_coords, lot_column, strain_column):
    """
    Find the lookup table for SR references

    :param worksheet: The worksheet object
    :params titer_coords: The coordinates of the titer block, might only need the last row
    :params lot_column: The column index for sr_lot
    :params strain_column: The column index for strain name
    :return: dictionary of mapping of serum ID to serum information for strain and lot ID

    """
    lot_mapping={}
    strain_mapping={}

    store=False
    for i in range(titer_coords['row_end']+1, worksheet.nrows, 1):
        if(worksheet.cell_value(i, 0) == "REFERENCE ANTISERUM"):
            print(f"Found REFERENCE ANTISERUM at row: {i}; col: 0")
            store=True

        if store and worksheet.cell_value(i, 0) !="":
            lot_mapping[worksheet.cell_value(i, 0)] = worksheet.cell_value(i, lot_column)
            if worksheet.cell_value(i, strain_column) !="":
                strain_mapping[worksheet.cell_value(i, 0)] = worksheet.cell_value(i, strain_column)
            else:
                if worksheet.cell_value(i, strain_column-1) !="":
                    strain_mapping[worksheet.cell_value(i, 0)] = worksheet.cell_value(i, strain_column-1)


    return {
        "lot_mapping": lot_mapping,
        "strain_mapping": strain_mapping
    }

def find_assay_date(worksheet, default="YYYY-MM-DD"):
    """
    Find the assay date

    :param worksheet: The worksheet object
    :return: the assay date
    """
    iso_date=default
    for i in range(worksheet.nrows):
        for j in range(worksheet.ncols):
            cell_value = str(worksheet.cell_value(i, j)).strip()
            if cell_value.startswith("DATE TESTED:"):
                date_str = cell_value.split(":", 1)[1].strip()
                date_obj = datetime.strptime(date_str, "%m/%d/%Y")
                iso_date = date_obj.date().isoformat()
                return iso_date

    return iso_date

def main():
    args = parse_args()

    # Pattern matching
    virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    virus_passage_pattern = r"(S\d+|E\d+|Spf|QMC)"
    serum_id_pattern = r"^[A-Z]$"
    serum_passage_pattern = r"(SIAT|EGG|SPA|QMC|CELL)"
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
            'col_end': titer_block["col_end"][0][0] + int(args.col_end_adjustment), # ADJUSTMENT! Specific to dataset, due two numeric columns to the right (4, titer value)
            'row_start': titer_block["row_start"][0][0],
            'row_end': titer_block["row_end"][0][0] # ADJUSTMENT! Specific to dataset
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
        lookups = find_sr_lookup(
            worksheet=worksheet,
            titer_coords=titer_coords,
            lot_column=int(args.lot_column),
            strain_column=virus_block['virus_col_idx']
            )

        lot_map=lookups["lot_mapping"]
        test_date=find_assay_date(worksheet=worksheet, default=args.test_date)

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
        # serum_mapping = serum_block["serum_mapping"]
        strain_map=lookups["strain_mapping"]

        print(f"Writing to file: {filebasename.replace(' ', '_')}_{worksheet.name.replace(' ', '_')}.tsv")

        test_subtype="None"
        with open(f"{filebasename.replace(' ', '_')}_{worksheet.name.replace(' ', '_')}.tsv", 'w') as outfile:
          # Print header to flattened titer values output file
          outfile.write("%s\n" % ("\t".join(["test_date","test_protocol","test_subtype","ag_strain_name","ag_passage","ag_cdc_id","sr_strain_name","sr_passage","sr_lot","titer_value"])))
          #outfile.write("%s\n" % ("\t".join(["test_date","assay_type","virus_strain","virus_strain_passage","serum_strain","serum_antigen_passage","serum_id","titer_value"])))
          for i in range(titer_coords['row_start'], titer_coords['row_end']):

              # Pull ag_strain information since it is consistent across a row.
              ag_strain_name=str(worksheet.cell_value(i,virus_block['virus_col_idx'])).strip().replace(" (NEW)","").replace("(NEW)", "").replace(" NEW","")
              ag_passage=str(worksheet.cell_value(i,virus_block['virus_passage_col_idx'])).strip().split("(", 1)[0].strip()
              ag_cdc_id=""

              #ag_cdc_id_col=virus_block['virus_col_idx']+1
              ag_cdc_id_col=int(args.cdcid_column)
              cell = worksheet.cell(i, ag_cdc_id_col)
              if cell.ctype == xlrd.XL_CELL_NUMBER:
                  ag_cdc_id = str(int(cell.value)) if cell.value.is_integer() else str(cell.value)
              else:
                  ag_cdc_id = str(cell.value).strip()

              if str(worksheet.cell_value(i,virus_block['virus_col_idx']+4)).strip()!="":
                  subtype_col=virus_block['virus_col_idx']+4
                  subtype_col=int(args.subtype_column)
                  test_subtype=str(worksheet.cell_value(i,subtype_col)).strip().replace("H3N2","H3")

              for j in range(titer_coords['col_start'], titer_coords['col_end']):
                  if not is_titer_value(worksheet.cell_value(i,j)):
                      break

                  # Pull sr_strain information for the particular cell in a row
                  #sr_strain_abbr=str(worksheet.cell_value(serum_block['serum_abbrev_row_idx'],j)).strip()
                  #sr_strain_name=serum_mapping.get(sr_strain_abbr, sr_strain_abbr)
                  # Turns out human pool is only listed in the REFERENCE ANTISERUM section
                  sr_strain_name=strain_map[str(worksheet.cell_value(serum_block['serum_id_row_idx'],j)).strip()].replace(" (NEW)","").replace("(NEW)", "").replace(" NEW","")
                  sr_passage=str(worksheet.cell_value(serum_block['serum_passage_row_idx'],j)).strip()
                  sr_lot=lot_map[str(worksheet.cell_value(serum_block['serum_id_row_idx'],j)).strip()]

                  if sr_strain_name=="Human POOL / CSID 3100021731,3100021827,3100021837":
                      sr_strain_name="A/Darwin/6/2021"
                      sr_passage="CELL"

                  # Pull titer value
                  titer = str(worksheet.cell_value(i,j)).strip()

                  # Print flattened titer values to output file
                  outfile.write(f"{test_date}\t{args.test_protocol}\t{test_subtype}\t{ag_strain_name}\t{ag_passage}\t{ag_cdc_id}\t{sr_strain_name}\t{sr_passage}\t{sr_lot}\t{titer}\n")

if __name__ == "__main__":
    main()
