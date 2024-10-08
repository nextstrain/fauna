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
from typing import Iterator, Optional, Tuple
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.flu_upload import flu_upload
from titer_block import find_titer_block, find_serum_rows, find_virus_columns

parser.add_argument('--assay_type', default='hi')
parser.add_argument('--human-ref-only', action="store_true",
    help="Only ingest human sera references, used for backfilling data that was skipped in previous ingests.")

ELIFE_COLUMNS = ["virus_strain", "serum_strain","serum_id", "titer", "source", "virus_passage", "virus_passage_category", "serum_passage", "serum_passage_category", "assay_type"]
EXPECTED_SUBTYPES = {"h1n1pdm", "h3n2", "vic", "yam"}
HUMAN_SERA_YEAR_REGEX = r"SH(vax|VAX|\s)?(\d{4})"

# Vaccine mapping used for mapping human pooled sera to a specific reference virus
# This is based on the proxy reference viruses used for human pooled sera
# in the flat files provided by VIDRL.
VACCINE_MAPPING = {
    "2023": {
        "egg": {
            "h1n1pdm": "A/Sydney/5/2021",
            "h3n2": "A/Darwin/9/2021",
            "vic": "B/Austria/1359417/2021",
            "yam": "B/Phuket/3073/2013"
        },
        "cell": {
            "h1n1pdm": "A/Sydney/5/2021",
            "h3n2": "A/Darwin/6/2021",
            "vic": "B/Austria/1359417/2021",
            "yam": "B/Phuket/3073/2013"
        }
    },
    "2024": {
        "egg": {
            "h1n1pdm": "A/Victoria/4897/2022",
            "h3n2": "A/Thailand/8/2022",
            "vic": "B/Austria/1359417/2021",
            "yam": "B/Phuket/3073/2013"
        },
        "cell": {
            "h1n1pdm": "A/Victoria/4897/2022",
            "h3n2": "A/Thailand/8/2022",
            "vic": "B/Austria/1359417/2021",
            "yam": "B/Phuket/3073/2013"
        }
    }
}

def parse_human_serum_id(original_id: str, year_regex: str) -> Tuple[Optional[str],Optional[str]]:
    """
    Attempts to parse the year (YYYY) from the provided *original_id* to
    construct a new standard serum id with the provided *year_regex*.

    Returns None if the year could not be parsed from the *original_id*.
    If year is successfully parsed, then returns a tuple of (year, new_serum_id).
    """
    year = new_serum_id = None
    matches = re.match(year_regex, original_id)
    if matches is None:
        return (None, None)

    year = matches.group(2)
    # Follow a standard pattern where serum_id is `Human pool <year>`
    # Need "human" in serum_id because this is how we match for human sera in seasonal flu
    # <https://github.com/nextstrain/seasonal-flu/blob/89f6cfd11481b2c51c50d68822c18d46ed56db51/workflow/snakemake_rules/download_from_fauna.smk#L93>
    new_serum_id = f"Human pool {year}"
    return (year, new_serum_id)


def parse_human_serum_references(human_serum_data, subtype):
    """
    Expects the *human_serum_data* from titer_block.find_serum_rows
    Returns parsed human serum references, where keys are the column number of
    the human serum reference in the Excel sheet and the values are the serum
    data with serum id, serum passage, and serum strain.
    """
    human_serum_references = {}
    egg_or_cell_regex = r"^(egg|cell)$" # Used with re.IGNORECASE

    potential_year_fields = ['serum_id', 'serum_passage', 'serum_abbrev']
    potential_egg_or_cell_fields = ['serum_passage', 'extra_info']

    for human_serum in human_serum_data:
        column = human_serum['col_idx']
        # First try to parse the year from the human serum data
        year = new_serum_id = None
        for field in potential_year_fields:
            year, new_serum_id = parse_human_serum_id(human_serum[field], HUMAN_SERA_YEAR_REGEX)
            # Use the first match of the potential fields
            if new_serum_id is not None:
                break

        # year is required to know which vaccine reference strain to use
        # Raise an error because this info should _always_ be available
        if year is None:
            raise Exception(f"Unable to process human sera column {column} ",
                            f"because none of {potential_year_fields} fields ",
                            f"matched the year regex {HUMAN_SERA_YEAR_REGEX!r}")

        # Then try to parse egg or cell from the human serum data
        egg_or_cell = None
        for field in potential_egg_or_cell_fields:
            matches = re.match(egg_or_cell_regex, human_serum[field], re.IGNORECASE)
            # Use the first match of the potential fields
            if matches is not None:
                egg_or_cell = matches.group(1).lower()
                break

        # egg_or_cell is required to know which vaccine reference strain to use,
        # so skip the human serum if it can't be parsed
        # Only outputting a warning because I've seen Excel worksheets _without_
        # any egg/cell distinctions from 2023. This will require extra correspondence
        # with VIDRL, so don't let it block ingest of other data.
        #   -Jover, 28 August 2024
        if egg_or_cell is None:
            print(f"WARNING: Skipping human sera column {column} ",
                  f"because none of {potential_egg_or_cell_fields} fields ",
                  f"matched the regex {egg_or_cell_regex!r}")
            continue

        # Raise a loud error so we know to update the VACCINE_MAPPING as needed
        try:
            serum_strain = VACCINE_MAPPING[year][egg_or_cell][subtype]
        except KeyError as err:
            raise Exception(f"VACCINE_MAPPING needs to be updated!") from err

        human_serum_references[column] = {
            "serum_id": new_serum_id,
            "serum_passage": egg_or_cell,
            "serum_strain": serum_strain
        }

    return human_serum_references


def read_vidrl(path, fstem, assay_type, subtype, human_ref_only):
    '''
    Read all csv tables in path, create data frame with reference viruses as columns
    '''
    exten = [ os.path.isfile(path + fstem + ext) for ext in ['.xls', '.xlsm', '.xlsx'] ]

    if True in exten:
        ind = exten.index(True)
        convert_vidrl_xls_to_tsv(path, fstem, ind, assay_type, subtype, human_ref_only)
    else:
        print("Unable to recognize file {}/{}".format(path,fstem))
        sys.exit()

def convert_vidrl_xls_to_tsv(path, fstem, ind, assay_type, subtype, human_ref_only):
    exts = ['.xls', '.xlsm', '.xlsx']
    workbook = xlrd.open_workbook(path+fstem + exts[ind])

    # Default patterns, VIDRL
    virus_pattern = r"[A-Z]/[\w\s-]+/.+/\d{4}"
    virus_passage_pattern = r"(MDCK|SIAT|E\d+|hCK)"
    serum_id_pattern = r"^[A-Z]\d{4,8}"
    serum_passage_pattern = r"(MDCK\d+|SIAT\d+|E\d+)"
    serum_abbrev_pattern = r"\w+\s{0,1}\w+/\d+.*"
    human_serum_pattern = r"(^SH\d+|SHVAX|SHvax|sera|vaxpool).*"
    crick = False

    for worksheet_index, worksheet in enumerate(workbook.sheets(), start=1):
        print(f"Reading worksheet {worksheet_index} '{worksheet.name}' in file '{fstem}'")
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
            ignore_serum_pattern=human_serum_pattern,
            log_human_sera=True,
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
        # # Use manual correction when needed
        # serum_mapping = {
        #     'Dar/6': 'A/Darwin/6/2021',
        #     'Thai/8': 'A/Thailand/8/2022',
        #     'Syd/856': 'A/Sydney/856/2023',
        #     'Can/327': 'A/Canberra/327/2023',
        #     'Can/302': 'A/Canberra/302/2023',
        #     'Vic/2997': 'A/Victoria/2997/2023',
        #     'Can/373': 'A/Canberra/373/2023',
        #     'Can/309': 'A/Canberra/309/2023',
        # }
        # print(f"corrected: serum_mapping={json.dumps(serum_mapping, indent=4)}")

        human_serum_references = parse_human_serum_references(serum_block['human_serum_data'], args.subtype)

        print("Human pooled serum references parsed from serum block")
        for col, values in human_serum_references.items():
            print(f"Column {col!r}: {values}")

        # Check if all the necessary indices were found
        if virus_block["virus_col_idx"] is None:
            print(f"Virus column index not found. Check the virus pattern: '{virus_pattern}'")
            break

        if virus_block["virus_passage_col_idx"] is None:
            print(f"Virus passage column index not found. Check the virus passage pattern: '{virus_passage_pattern}'")
            break

        if serum_block["serum_id_row_idx"] is None:
            print(f"Serum ID row index not found. Check the serum ID pattern: '{serum_id_pattern}'")
            break

        if serum_block["serum_passage_row_idx"] is None:
            print(f"Serum passage row index not found. Check the serum passage pattern: '{serum_passage_pattern}'")
            break

        if serum_block["serum_abbrev_row_idx"] is None:
            print(f"Serum abbreviated name row index not found. Check the serum abbreviated name pattern: '{serum_abbrev_pattern}'")
            break

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
                    # Special handling of human pooled sera that were matched to
                    # vaccine reference strain instead of the normal serum mapping
                    if j in human_serum_references:
                        serum_id = human_serum_references[j]['serum_id']
                        serum_passage = human_serum_references[j]['serum_passage']
                        serum_strain = human_serum_references[j]['serum_strain']
                    # Skip other titer measurements if we only want to ingest human serum references
                    elif human_ref_only:
                        continue
                    else:
                        serum_id = str(mat.cell_value(serum_id_row_index,j)).strip().replace(' ','')
                        serum_passage = str(mat.cell_value(serum_passage_row_index,j)).strip()
                        serum_abbr = str(mat.cell_value(serum_strain_row_index,j)).strip()
                        serum_abbr = serum_abbr.replace(' ','')
                        serum_strain = serum_mapping.get(serum_abbr, serum_abbr)

                    titer = str(mat.cell_value(i,j)).strip()
                    line = "%s\n" % ("\t".join([ virus_strain, serum_strain, serum_id, titer, source, virus_passage, virus_passage_category, serum_passage, serum_passage_category, assay_type]))
                    outfile.write(line)


def read_csv_to_dict(csv_file: str) -> Iterator[dict]:
    """
    Read provided *csv_file* and yields each row as a dict.
    """
    with open(csv_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=",")
        yield from reader


def curate_flat_records(records: Iterator[dict], fstem: str, assay_type: str) -> Iterator[dict]:
    """
    Curate the measurement records expected from the _flat_file.csv files.
    """
    # The new column names need to be one of the ELIFE_COLUMNS in order to be
    # included in the temporary output file that's then passed to elife_upload.py
    column_map = {
        "test virus": "virus_strain",
        "test virus passage": "virus_passage",
        "reference antigen": "serum_strain",
        "antisera passage": "serum_passage",
        "ferret": "serum_id",
        "titre": "titer",
        "antisera": "serum_abbr",
    }
    for record in records:
        new_record = {new_field: record[old_field] for old_field, new_field in column_map.items()}
        new_record["assay_type"] = assay_type
        new_record["virus_passage_category"] = ""
        new_record["serum_passage_category"] = ""
        new_record["source"] = "vidrl_{}.csv".format(fstem)

        # Standardized human sera name and ids
        if new_record["serum_id"] == "NA":
            human_serum_id = new_record["serum_abbr"]
            egg_or_cell = human_serum_id[-1].lower()

            if egg_or_cell == "e":
                new_record["serum_passage"] = "egg"
            elif egg_or_cell == "c":
                new_record["serum_passage"] = "cell"
            else:
                raise ValueError(f"Human serum id {human_serum_id!r} does not include expected egg or cell type")

            # TODO: Verify the serum strain matches the year's vaccine strain
            new_record["serum_strain"] = re.sub(r"pool$", "", new_record["serum_strain"])
            _, new_record["serum_id"] = parse_human_serum_id(human_serum_id, HUMAN_SERA_YEAR_REGEX)

        yield new_record


def write_records_to_tsv(records: Iterator[dict], output_file: str):
    with open(output_file, "w", newline="") as fh:
        tsv_writer = csv.DictWriter(
            fh,
            ELIFE_COLUMNS,
            extrasaction="ignore",
            delimiter="\t",
            lineterminator="\n",
        )
        tsv_writer.writeheader()
        for record in records:
            tsv_writer.writerow(record)


def read_flat_vidrl(path, fstem, assay_type):
    """
    Read the flat CSV file with *fstem* in the provided *path* and convert
    to the expected TSV file at `data/tmp/<fstem>.tsv` for tdb/elife_upload.
    """

    filepath = path + fstem + ".csv"
    output_filepath ="data/tmp/{}.tsv".format(fstem)

    records = read_csv_to_dict(filepath)
    curated_records = curate_flat_records(records, fstem, assay_type)
    write_records_to_tsv(curated_records, output_filepath)


if __name__=="__main__":
    args = parser.parse_args()
    # Asserting here because this is using a shared parser
    # other tdb scripts do not require subtype
    assert args.subtype is not None, "Subtype needs to be specified with --subtype"
    assert args.subtype in EXPECTED_SUBTYPES, f"Subtype must be one of {EXPECTED_SUBTYPES!r}"

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
        read_vidrl(args.path, args.fstem, args.assay_type, args.subtype, args.human_ref_only)

    if args.preview:
        command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + args.subtype + " --path data/tmp/ --fstem " + args.fstem + " --preview"
        print(command)
        subprocess.call(command, shell=True)
    else:
        command = "python tdb/elife_upload.py -db " + args.database +  " --subtype " + args.subtype + " --path data/tmp/ --fstem " + args.fstem
        print(command)
        subprocess.call(command, shell=True)
