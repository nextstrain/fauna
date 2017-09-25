from __future__ import print_function
import os, sys
from pdb import set_trace


"""
This script takes a TSV file with the following fields:
0: accession
1: sequence name
2: date (in nextstrain 2016-01-20 format)
3: host
4: country
5: division (admin1, i.e. region / provence / state)
6: genotype
7: whole genome sequence

And turns it into a ViPR-like fasta file for uploding through the fauna architecture

mumps vipr-like format:
0: accession, 1: strain, 2: date, 3: host, 4: country, 5: division (admin1), 6: muv_genotype
"""

if __name__=="__main__":
    with open(sys.argv[1], 'rU') as fh:
        for line in fh:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            try:
                assert(len(fields) == 8)
            except AssertionError:
                print("Line {} didn't have 8 fields!".format(line))
            print(">" + '|'.join(fields[0:7]))
            print(fields[7])
