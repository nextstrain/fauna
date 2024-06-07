#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import xlrd
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description="Find the block of titers in an Excel worksheet.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--file",
                        default="~/nextstrain/fludata/VIDRL-Melbourne-WHO-CC/raw-data/A/H1N1pdm/HI/2024/20240528H1N1.xlsx",
                        required=False,
                        help="Path to the Excel file")
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
        if value.startswith('<'):
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
            next_cell_value = worksheet.cell_value(row_idx, col_idx + 1) if col_idx + 1 < worksheet.ncols else None

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
            next_cell_value = worksheet.cell_value(row_idx + 1, col_idx) if row_idx + 1 < worksheet.nrows else None

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
        "row_end": sorted_row_end
    }

def main():
    args = parse_args()

    # Load the Excel file
    workbook = xlrd.open_workbook(args.file)
    worksheet = workbook.sheet_by_index(0) # first sheet or loop

    # Find the block of titers in the worksheet
    titer_block = find_titer_block(worksheet)
    
    # Print the most likely row and column indices for the titer block
    print(f"Most likely (n={titer_block['col_start'][0][1]}) col_start: {titer_block['col_start'][0][0]}")
    print(f"Most likely (n={titer_block['col_end'][0][1]}) col_end: {titer_block['col_end'][0][0]}")
    print(f"Most likely (n={titer_block['row_start'][0][1]}) row_start: {titer_block['row_start'][0][0]}")
    print(f"Most likely (n={titer_block['row_end'][0][1]}) row_end: {titer_block['row_end'][0][0]}")

    ## TODO: Add logic to parse serum and virus annotations row and column indices
    ## TODO: Match abbreviated names across the top to the full names along the left side and auto convert to full namess

if __name__ == "__main__":
    main()
