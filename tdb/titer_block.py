#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import xlrd
from collections import defaultdict
import re


def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Find the block of titers in an Excel worksheet.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--file",
        default="~/nextstrain/fludata/VIDRL-Melbourne-WHO-CC/raw-data/A/H1N1pdm/HI/2024/20240528H1N1.xlsx",
        required=False,
        help="Path to the Excel file",
    )
    parser.add_argument(
        "--source",
        default="vidrl",
        required=False,
        help="Either vidrl, niid, or crick. Used to select appropriate pattern"
    )
    return parser.parse_args()

def is_numeric(value):
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


def find_titer_block(worksheet):
    """
    Find the block of titers in the worksheet.

    To reduce false positives, the function looks for rows and columns where at least two consecutive cells contain numeric values.

    :param worksheet: The worksheet object
    :return: Dictionary containing the most likely row and column indices for the titer block sorted by vote count
    """
    col_start_dict = defaultdict(int)
    col_end_dict = defaultdict(int)
    row_start_dict = defaultdict(int)
    row_end_dict = defaultdict(int)

    # Iterate over each row
    for row_idx in range(worksheet.nrows):
        first_numeric_index = None
        last_numeric_index = None
        for col_idx in range(worksheet.ncols):
            cell_value = worksheet.cell_value(row_idx, col_idx)
            next_cell_value = (
                worksheet.cell_value(row_idx, col_idx + 1)
                if col_idx + 1 < worksheet.ncols
                else None
            )

            # Check for two consecutive numeric values in a row before incrementing the count
            if is_numeric(cell_value) and is_numeric(next_cell_value):
                if first_numeric_index is None:
                    first_numeric_index = col_idx
                last_numeric_index = col_idx + 1

        # Only increment if any numeric values were found in the row
        if first_numeric_index is not None:
            col_start_dict[first_numeric_index] += 1
            col_end_dict[last_numeric_index] += 1

    # Iterate over each column
    for col_idx in range(worksheet.ncols):
        first_numeric_index = None
        last_numeric_index = None
        for row_idx in range(worksheet.nrows):
            cell_value = worksheet.cell_value(row_idx, col_idx)
            next_cell_value = (
                worksheet.cell_value(row_idx + 1, col_idx)
                if row_idx + 1 < worksheet.nrows
                else None
            )

            # Check for two consecutive numeric values in a column before incrementing the count
            if is_numeric(cell_value) and is_numeric(next_cell_value):
                if first_numeric_index is None:
                    first_numeric_index = row_idx
                last_numeric_index = row_idx + 1

        # Only increment if any numeric values were found in the column
        if first_numeric_index is not None:
            row_start_dict[first_numeric_index] += 1
            row_end_dict[last_numeric_index] += 1

    # Sort the dictionaries by frequency in descending order
    sorted_col_start = sorted(col_start_dict.items(), key=lambda item: item[1], reverse=True)
    sorted_col_end = sorted(col_end_dict.items(), key=lambda item: item[1], reverse=True)
    sorted_row_start = sorted(row_start_dict.items(), key=lambda item: item[1], reverse=True)
    sorted_row_end = sorted(row_end_dict.items(), key=lambda item: item[1], reverse=True)

    return {
        "col_start": sorted_col_start,
        "col_end": sorted_col_end,
        "row_start": sorted_row_start,
        "row_end": sorted_row_end,
    }


def find_virus_columns(worksheet, titer_coords, virus_pattern=None, virus_passage_pattern=None):
    """
    Find the columns containing virus names based on the most likely column indices for the titer block.

    :param worksheet: The worksheet object
    :param titer_block: Dictionary containing the titer block coordinates
    :param virus_pattern: Regular expression pattern to match virus names (optional)
    :param virus_passage_pattern: Regular expression pattern to match virus passage data (optional)
    :return: Dictionary containing the column indices for virus names, virus passage data, and a list of virus names
    """
    # List of outputs
    virus_col_idx = None  # Index of the column containing virus names
    virus_passage_col_idx = None  # Index of the column containing virus passage data
    virus_names = (
        []
    )  # List of virus names, will be used by find_antigen_rows to map abbreviated antigen names to full names

    # By default use VIDRL-Melbourne-WHO-CC specific patterns
    # Define a regular expression pattern to match virus names
    if virus_pattern is None:
        virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    # Define a regular expression pattern to match virus passage column
    if virus_passage_pattern is None:
        virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"

    # Find the column containing virus names searching to the left of the titer block
    for col_idx in range(titer_coords['col_start'] - 1, -1, -1):
        virus_count = 0
        for row_idx in range(titer_coords['row_start'], titer_coords['row_end'] + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.search(virus_pattern, cell_value):
                virus_count += 1

        # Index of the first column that contains more than 50% rows matching the virus pattern
        if virus_count > (titer_coords['row_end'] - titer_coords['row_start']) / 2:
            virus_col_idx = col_idx
            break

    # Get the virus names from the column containing virus names
    # This allows for some lienency in matching the virus pattern in the column
    for row_idx in range(titer_coords['row_start'], titer_coords['row_end'] + 1):
        cell_value = str(worksheet.cell_value(row_idx, virus_col_idx))
        if r"/" in cell_value:
            virus_names.append(cell_value)

    # Find the column containing virus passage data searching to the right of the titer block
    # If not found, search left
    for col_idx in list(range(titer_coords['col_end'], worksheet.ncols)) + list(range(titer_coords['col_start'] - 1, -1, -1)):
        virus_passage_count = 0
        for row_idx in range(titer_coords['row_end']):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.search(virus_passage_pattern, cell_value):
                virus_passage_count += 1

        # Index of the first column that contains more than 50% rows matching the virus passage pattern
        if virus_passage_count > (titer_coords['row_end'] - titer_coords['row_start']) / 2:
            virus_passage_col_idx = col_idx
            break

    return {
        "virus_col_idx": virus_col_idx,
        "virus_passage_col_idx": virus_passage_col_idx,
        "virus_names": virus_names,
    }


def find_serum_rows(worksheet, titer_coords, virus_names=None, serum_id_pattern=None, serum_passage_pattern=None, serum_abbrev_pattern=None, crick=False, ignore_serum_pattern=None):
    """
    Find the row containing cell passage data and the row containing abbreviated serum names.

    :param worksheet: The worksheet object
    :param titer_coords: Dictionary containing the titer block coordinates
    :param virus_names: List of virus names
    :param serum_id_pattern: Regular expression pattern to match serum ID (optional)
    :param serum_passage_pattern: Regular expression pattern to match cell passage data (optional)
    :param serum_abbrev_pattern: Regular expression pattern to match abbreviated serum names (optional)
    :return: Dictionary containing the row indices for serum ID, cell passage data, abbreviated serum names, and a mapping of abbreviated names to full names
    """
    # List of outputs
    serum_id_row_idx = None  # Index of the row containing serum ID
    serum_passage_row_idx = None  # Index of the row containing cell passage data
    serum_abbrev_row_idx = None  # Index of the row containing abbreviated antigen names
    serum_mapping = {}  # Mapping of abbreviated antigen names to full names

    # By default, use VIDRL-Melbourne-WHO-CC specific patterns
    # Define a regular expression pattern to match serum ID
    if serum_id_pattern is None:
        serum_id_pattern = r"^[A-Z]\d{4,8}$"
    # Define a regular expression pattern to match cell passage data
    if serum_passage_pattern is None:
        serum_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+)"
    # Define a regular expression pattern to match abbreviated antigen names
    if serum_abbrev_pattern is None:
        serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"
    if ignore_serum_pattern is None:
        ignore_serum_pattern = r"(^SH\d+|SHVAX|SHvax|sera).*"

    # Find the row containing serum ID searching from the top of the titer block upwards
    for row_idx in range(titer_coords['row_start'] - 1, -1, -1):  # Iterate from row_start to the top
        serum_id_count = 0
        for col_idx in range(titer_coords['col_start'], titer_coords['col_end'] + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            cell_value = re.sub(r'[\r\n ]+', '', cell_value)
            if re.search(serum_id_pattern, cell_value):
                serum_id_count += 1
        # Index of the first row that contains more than 50% columns matching the serum ID pattern
        if serum_id_count > (titer_coords['col_end'] - titer_coords['col_start']) / 2:
            serum_id_row_idx = row_idx
            break


    # Find the row containing cell passage data searching from the top of the titer block upwards
    for row_idx in range(titer_coords['row_start'] - 1, -1, -1):  # Iterate from row_start to the top
        serum_passage_count = 0
        for col_idx in range(titer_coords['col_start'], titer_coords['col_end'] + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            cell_value = re.sub(r'[\r\n ]+', '', cell_value)
            if re.search(serum_passage_pattern, cell_value):
                serum_passage_count += 1
        if serum_passage_count > (titer_coords['col_end'] - titer_coords['col_start']) / 2:
            serum_passage_row_idx = row_idx
            break

    # Find the row containing abbreviated serum names searching from the top of the titer block upwards
    for row_idx in range(titer_coords['row_start'] -1, -1, -1):
        serum_abbrev_count = 0
        for col_idx in range(titer_coords['col_start'], titer_coords['col_end'] + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            cell_value = re.sub(r'[\r\n ]+', '', cell_value)
            if re.search(serum_abbrev_pattern, cell_value):
                serum_abbrev_count += 1

        if serum_abbrev_count > (titer_coords['col_end'] - titer_coords['col_start']) / 2:
            serum_abbrev_row_idx = row_idx
            break

    # Map abbreviated serum names to full names
    virus_idx = 0
    for col_idx in range(titer_coords['col_start'], titer_coords['col_end'] + 1):
        cell_value = str(worksheet.cell_value(serum_abbrev_row_idx, col_idx))
        if crick:
            cell_value = cell_value + str(worksheet.cell_value(serum_abbrev_row_idx+1, col_idx))
        cell_value = re.sub(r'[\r\n ]+', '', cell_value)
        if cell_value == "":
            break
        # Ignore human serum (e.g. "SH2002", "sera", "SHVAX2002")
        if re.search(ignore_serum_pattern, cell_value):
            break
            serum_mapping[cell_value] = virus_names[virus_idx]
            virus_idx += 1

    return {
        "serum_id_row_idx": serum_id_row_idx,
        "serum_passage_row_idx": serum_passage_row_idx,
        "serum_abbrev_row_idx": serum_abbrev_row_idx,
        "serum_mapping": serum_mapping,
    }


def main():
    args = parse_args()

    # Default patterns, VIDRL
    virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"
    serum_id_pattern = r"^[A-Z]\d{4,8}$"
    serum_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+)"
    serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"
    crick = False

    if(args.source == "niid"):
        serum_id_pattern = r".+(No\.|no\.).+"
        serum_passage_pattern = r".+(Egg|Cell).+"

    if(args.source == "crick"):
        virus_pattern = r"[A-Z]/[\w\s-]+"
        serum_id_pattern = r"F\d+/\d+"
        serum_passage_pattern = r"(MDCK|SIAT|Egg)"
        serum_abbrev_pattern = r"[A-Z]/[\w\s-]+"
        crick = True

    # Load the Excel file
    workbook = xlrd.open_workbook(args.file)

    for worksheet in workbook.sheets():
        print(f"Worksheet: {worksheet.name}")

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
            # print("Serum mapping:")
            # for abbrev, full in serum_block["serum_mapping"].items():
            #     print(f"  {abbrev} -> {full}")

            print("serum_mapping = {")
            for abbrev, full in serum_block["serum_mapping"].items():
                print(f"    '{abbrev}': '{full}',")
            print("}")


if __name__ == "__main__":
    main()
