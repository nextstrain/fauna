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
    return parser.parse_args()


def is_numeric(value):
    """
    Check if the value is numeric or a string representing a numeric value with '<'.

    Basically checking for titer values e.g. "80", "< 10", "1400", etc.
    """
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, str):
        value = value.strip()
        if value.startswith("<"):
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


def find_virus_columns(worksheet, col_start, col_end, row_start, row_end):
    """
    Find the columns containing virus names based on the most likely column indices for the titer block.
    """
    # List of outputs
    virus_col_idx = None  # Index of the column containing virus names
    virus_passage_col_idx = None  # Index of the column containing virus passage data
    virus_names = (
        []
    )  # List of virus names, will be used by find_antigen_rows to map abbreviated antigen names to full names

    # VIDRL-Melbourne-WHO-CC specific patterns (TODO: parameterize the pattern matching by data source)
    # Define a regular expression pattern to match virus names
    virus_pattern = r"[A-Z]/\w+/.+/\d{4}"
    # Define a regular expression pattern to match virus passage column
    virus_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+)"

    # Find the column containing virus names searching to the left of the titer block
    for col_idx in range(col_start - 1, -1, -1):
        virus_count = 0
        for row_idx in range(row_start, row_end + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.match(virus_pattern, cell_value):
                virus_count += 1

        # Index of the first column that contains more than 50% rows matching the virus pattern
        if virus_count > (row_end - row_start) / 2:
            virus_col_idx = col_idx
            break

    # Get the virus names from the column containing virus names
    # This allows for some lienency in matching the virus pattern in the column
    for row_idx in range(row_start, row_end + 1):
        virus_names.append(str(worksheet.cell_value(row_idx, virus_col_idx)))

    # Find the column containing virus passage data searching to the right of the titer block
    for col_idx in range(col_end, worksheet.ncols):
        virus_passage_count = 0
        for row_idx in range(row_end):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.match(virus_passage_pattern, cell_value):
                virus_passage_count += 1

        # Index of the first column that contains more than 50% rows matching the virus passage pattern
        if virus_passage_count > (row_end - row_start) / 2:
            virus_passage_col_idx = col_idx
            break

    return {
        "virus_col_idx": virus_col_idx,
        "virus_passage_col_idx": virus_passage_col_idx,
        "virus_names": virus_names,
    }


def find_serum_rows(worksheet, col_start, col_end, row_start, row_end, virus_names=None):
    """
    Find the row containing cell passage data and the row containing abbreviated serum names.
    """
    # List of outputs
    serum_id_row_idx = None  # Index of the row containing serum ID
    serum_passage_row_idx = None  # Index of the row containing cell passage data
    serum_abbrev_row_idx = None  # Index of the row containing abbreviated antigen names
    serum_mapping = {}  # Mapping of abbreviated antigen names to full names

    # VIDRL-Melbourne-WHO-CC specific patterns (TODO: parameterize the pattern matching by data source)
    # Define a regular expression pattern to match serum ID
    serum_id_pattern = r"^[A-Z]\d{4,8}$"
    # Define a regular expression pattern to match cell passage data
    serum_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+)"
    # Define a regular expression pattern to match abbreviated antigen names
    serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"

    # Find the row containing serum ID searching from the top of the titer block upwards
    for row_idx in range(row_start - 1, -1, -1):  # Iterate from row_start to the top
        serum_id_count = 0
        for col_idx in range(col_start, col_end + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.match(serum_id_pattern, cell_value):
                serum_id_count += 1
        # Index of the first row that contains more than 50% columns matching the serum ID pattern
        if serum_id_count > (col_end - col_start) / 2:
            serum_id_row_idx = row_idx
            break

    # Find the row containing cell passage data searching from the top of the titer block upwards
    for row_idx in range(row_start - 1, -1, -1):  # Iterate from row_start to the top
        serum_passage_count = 0
        for col_idx in range(col_start, col_end + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.match(serum_passage_pattern, cell_value):
                serum_passage_count += 1
        if serum_passage_count > (col_end - col_start) / 2:
            serum_passage_row_idx = row_idx
            break

    # Find the row containing abbreviated serum names searching from the top of the titer block upwards
    for row_idx in range(row_start -1, -1, -1):
        serum_abbrev_count = 0
        for col_idx in range(col_start, col_end + 1):
            cell_value = str(worksheet.cell_value(row_idx, col_idx))
            if re.match(serum_abbrev_pattern, cell_value):
                serum_abbrev_count += 1

        if serum_abbrev_count > (col_end - col_start) / 2:
            serum_abbrev_row_idx = row_idx
            break

    # Map abbreviated serum names to full names
    virus_idx = 0
    for col_idx in range(col_start, col_end + 1):
        cell_value = str(worksheet.cell_value(serum_abbrev_row_idx, col_idx))
        # A more lenient check for the presence of a "/" in the cell value to find the abbreviated serum names
        if r"/" in cell_value:
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

    # Load the Excel file
    workbook = xlrd.open_workbook(args.file)
    worksheet = workbook.sheet_by_index(0)  # first sheet or loop

    # Find the block of titers in the worksheet
    titer_block = find_titer_block(worksheet)

    virus_block = find_virus_columns(
        worksheet=worksheet,
        col_start=titer_block["col_start"][0][0],
        col_end=titer_block["col_end"][0][0],
        row_start=titer_block["row_start"][0][0],
        row_end=titer_block["row_end"][0][0],
    )
    serum_block = find_serum_rows(
        worksheet=worksheet,
        col_start=titer_block["col_start"][0][0],
        col_end=titer_block["col_end"][0][0],
        row_start=titer_block["row_start"][0][0],
        row_end=titer_block["row_end"][0][0],
        virus_names=virus_block["virus_names"],
    )

    # Print the most likely row and column indices for the titer block
    print(f"Titer block: n = {titer_block['row_start'][0][1]}x{titer_block['col_start'][0][1]} = {titer_block['row_start'][0][1]*titer_block['col_start'][0][1]}")
    print(f"  Most likely (n={titer_block['col_start'][0][1]}) col_start: {titer_block['col_start'][0][0]}")
    print(f"  Most likely (n={titer_block['col_end'][0][1]}) col_end: {titer_block['col_end'][0][0]}")
    print(f"  Most likely (n={titer_block['row_start'][0][1]}) row_start: {titer_block['row_start'][0][0]}")
    print(f"  Most likely (n={titer_block['row_end'][0][1]}) row_end: {titer_block['row_end'][0][0]}")

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
