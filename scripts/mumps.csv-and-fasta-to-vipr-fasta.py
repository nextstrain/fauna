from __future__ import print_function
import os, sys
from pdb import set_trace
from Bio import SeqIO


"""
How to run:
python scripts/mumps.csv-and-fasta-to-vipr-fasta.py INPUT_FASTA INPUT_CSV OUTPUT_FASTA

Inputs:
(1) FASTA file with sequences

(2) CSV file with fields (no header, comment lines starting with #)
    1. name as found in FASTA input file
    2. name for output FASTA
    3. collection date (in nextstrain 2016-01-20 format)
    4. host
    5. country
    6. state/region
    7. MuV genotype

The output of this script is a FASTA file with headers displaying (seperated by '|')
    1. random numbering - this will later be filled in by GenBank accession
    2. strain name
    3. collection date
    4. host species
    5. country
    6. state/region
    7. genotype

P.S. if your CSV has ^M characters in it, run tr '\r' '\n' < mumps.csv > mumps.unix.csv

"""

if __name__=="__main__":
    # FASTA FILE
    records = list(SeqIO.parse(sys.argv[1], "fasta"))

    # CSV FILE
    db = {}
    count = 0
    with open(sys.argv[2], 'rU') as fh:
        for line in fh:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split(',')
            try:
                assert(len(fields) == 7)
            except AssertionError:
                print("Line {} didn't have 7 fields!".format(line))
                continue
            count += 1
            db[fields[0]] = "BCCDC{}|".format(count) + '|'.join(fields[1:7])

    # OUTPUT FASTA
    vipr = []
    for record in records:
        if record.name not in db:
            print("{} not in CSV. ignoring.".format(record.name))
            continue
        record.id = db[record.name]
        record.name = ''
        record.description = ''
        vipr.append(record)
    with open(sys.argv[3], "w") as fh:
        SeqIO.write(vipr, fh, "fasta")
