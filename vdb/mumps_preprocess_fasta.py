from __future__ import print_function
from Bio import SeqIO
import argparse
import os, sys, re
from pdb import set_trace

def fixHeaders(seqs):
    with open("source-data/mumps-name-fix.txt", 'rU') as fh:
        fixes = {x.strip().split()[0]: x.strip().split()[1] for x in fh.readlines() if not x.startswith('#')}
    fixed = []
    for seq in seqs:
        if seq.name in fixes:
            print("seq.name:", seq.name, " -> ", fixes[seq.name])
            seq.name = fixes[seq.name]
            fixed.append(seq)
        else:
            print("ignoring ", seq.name)
    print("Taking {} of {} sequences forward".format(len(fixed), len(seqs)))
    return fixed

def add_region(seqs):
    for seq in seqs:
        result = re.search('\|MuV[si]\/([A-z]+)\.',seq.name)
        if result:
            region = ''.join(result.groups())
        else:
            region = "NA"
        seq.name = seq.name + "|" + region
    return seqs

def standardiseViaName(seqs):
    for seq in seqs:
        seq.name = seq.name.replace('[','(').replace(']',')')
        seq.description = seq.name
        seq.id = seq.name
    return seqs

def vipr():
    with open("data/mumps_all_full.fasta", 'rU') as fh:
        seqs = [x for x in SeqIO.parse(fh, 'fasta')]
    seqs = fixHeaders(seqs)
    seqs = add_region(seqs)
    seqs = standardiseViaName(seqs)
    with open("data/mumps_all.fasta", 'w') as fh:
        SeqIO.write(seqs, fh, "fasta")
    print("Fixed VIPR fasta file saved to data/mumps_all.fasta")

def preprocessFASTA(fname):
    with open(fname, 'rU') as fh:
        seqs = [x for x in SeqIO.parse(fh, 'fasta')]
    seqs = fixHeaders(seqs)
    seqs = standardiseViaName(seqs)
    fname_out = fname + ".fixed.fasta"
    with open(fname_out, 'w') as fh:
        SeqIO.write(seqs, fh, "fasta")
    print("Fixed FASTA and saved to ", fname_out)

def collect_args():
    parser = argparse.ArgumentParser(description = "Preprocess mumps FASTA for fauna upload")
    parser.add_argument('--vipr', action='store_true', help='vipr fasta file')
    parser.add_argument('--fasta', type=str, default="", help='other fasta file')
    return parser.parse_args()

if __name__=="__main__":
    params = collect_args()
    if params.vipr:
        vipr()
    elif params.fasta:
        preprocessFASTA(params.fasta)

    # elif params.
