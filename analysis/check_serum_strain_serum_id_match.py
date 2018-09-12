# Open downloaded titer table and check matches between serum_strain and serum_id
# There should only be a single serum_strain for each serum_id
# There can be multiple serum_ids for each serum_strain

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('infile', default=None, type=str, help="file to test")

if __name__=="__main__":
    args = parser.parse_args()

    id_to_strain_mapping = {}

    if args.infile:
        with open(args.infile) as fd:
            rd = csv.reader(fd, delimiter="\t", quotechar='"')
            for row in rd:
                serum_strain = row[1]   # second row is serum serum_strain
                serum_id = row[2]       # third row is serum_id
                if serum_id in id_to_strain_mapping:
                    id_to_strain_mapping[serum_id].add(serum_strain)
                else:
                    id_to_strain_mapping[serum_id] = set([serum_strain])

    print("ALL SERUM_IDS")
    print(id_to_strain_mapping)
    print()

    print("PROBLEMATIC SERUM_IDS")
    for serum_id, serum_strains in id_to_strain_mapping.items():
        if len(serum_strains)>1:
            print("serum_id", serum_id)
            print("serum_strains", serum_strains)
