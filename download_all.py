# Simple script to run required operations to
# 1. Download FASTAs from database
# 2. Copy FASTAs to nextflu directory
# 3. Download titer tables from database
# 4. Copy titer tables to nextflu directory
# Run from base fauna directory with python flu/download_all.py
# Assumes that nextflu/, nextflu-cdc/ and nextflu-cdc-fra/ are
# sister directories to fauna/

import os, subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--copy_data', default=False, action="store_true", help="copy downloaded data to augur/data")
parser.add_argument('--sequences', default=False, action="store_true", help="download sequences from vdb")
parser.add_argument('--crick', default=False, action="store_true", help="download titers from crick reports")
parser.add_argument('--cdc', default=False, action="store_true", help="download titers from cdc reports")
parser.add_argument('--virus', default="flu", help="virus to download; default is flu")
parser.add_argument('--all', default=False, action="store_true", help="concatenate matched titer and strain TSVs")
parser.add_argument('--augur_path', default="../nextflu/augur/", help="path to the desired augur directory; default is ../nextflu/augur/")

if __name__=="__main__":
    args = parser.parse_args()

    if args.all:
        args.crick = True
        args.cdc = True

    if args.virus == "flu":
        # Download FASTAs from database
        if args.sequences:
            for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
                call = "python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_%s --fstem %s"%(lineage, lineage)
                os.system(call)

        # Copy FASTAs to nextflu/augur directory. Leave in fauna/data/ for nextstrain/augur.
        if args.copy_data:
            os.system("cp data/h3n2.fasta %sdata/"%(args.augur_path))
            os.system("cp data/h1n1pdm.fasta %sdata/"%(args.augur_path))
            os.system("cp data/vic.fasta %sdata/"%(args.augur_path))
            os.system("cp data/yam.fasta %sdata/"%(args.augur_path))

        if args.crick:
            for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
                call = "python tdb/download.py -db tdb -v flu --subtype %s --select assay_type:hi --fstem %s_crick_hi"%(lineage, lineage)
                os.system(call)


        # Download CDC HI titers from database
        if args.cdc:
            for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
                for passage in ["egg", "cell"]:
                    call = "python tdb/download.py -db cdc_tdb -v flu --subtype %s --select assay_type:hi serum_passage_category:%s --fstem %s_cdc_hi_%s"%(lineage, passage, lineage, passage)
                    os.system(call)

            # Download CDC FRA titers from database
            lineage = 'h3n2'
            for passage in ["egg", "cell"]:
                call = "python tdb/download.py -db cdc_tdb -v flu --subtype %s --select assay_type:fra serum_passage_category:%s --fstem %s_cdc_fra_%s"%(lineage, passage, lineage, passage)
                os.system(call)

        if args.all:

            # Concatenate strains TSVs for each subtype
            out = 'data/h3n2_all_hi_strains.tsv'
            hi_strains_h3n2 = ['data/h3n2_cdc_hi_cell_strains.tsv', 'data/h3n2_cdc_hi_egg_strains.tsv', 'data/h3n2_crick_hi_strains.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_strains_h3n2
                print call
                subprocess.call(call, stdout=f)

            out = 'data/h1n1pdm_all_hi_strains.tsv'
            hi_strains_h1n1pdm = ['data/h1n1pdm_cdc_hi_cell_strains.tsv', 'data/h1n1pdm_cdc_hi_egg_strains.tsv', 'data/h1n1pdm_crick_hi_strains.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_strains_h1n1pdm
                print call
                subprocess.call(call, stdout=f)

            out = 'data/vic_all_hi_strains.tsv'
            hi_strains_vic = ['data/vic_cdc_hi_cell_strains.tsv', 'data/vic_cdc_hi_egg_strains.tsv', 'data/vic_crick_hi_strains.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_strains_vic
                print call
                subprocess.call(call, stdout=f)

            out = 'data/yam_all_hi_strains.tsv'
            hi_strains_yam = ['data/yam_cdc_hi_cell_strains.tsv', 'data/yam_cdc_hi_egg_strains.tsv', 'data/yam_crick_hi_strains.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_strains_yam
                print call
                subprocess.call(call, stdout=f)


            # Concatenate titers TSVs for each subtype
            out = 'data/h3n2_all_hi_titers.tsv'
            hi_titers_h3n2 = ['data/h3n2_cdc_hi_cell_titers.tsv', 'data/h3n2_cdc_hi_egg_titers.tsv', 'data/h3n2_crick_hi_titers.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_titers_h3n2
                print call
                subprocess.call(call, stdout=f)

            out = 'data/h1n1pdm_all_hi_titers.tsv'
            hi_titers_h1n1pdm = ['data/h1n1pdm_cdc_hi_cell_titers.tsv', 'data/h1n1pdm_cdc_hi_egg_titers.tsv', 'data/h1n1pdm_crick_hi_titers.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_titers_h1n1pdm
                print call
                subprocess.call(call, stdout=f)

            out = 'data/vic_all_hi_titers.tsv'
            hi_titers_vic = ['data/vic_cdc_hi_cell_titers.tsv', 'data/vic_cdc_hi_egg_titers.tsv', 'data/vic_crick_hi_titers.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_titers_vic
                print call
                subprocess.call(call, stdout=f)

            out = 'data/yam_all_hi_titers.tsv'
            hi_titers_yam = ['data/yam_cdc_hi_cell_titers.tsv', 'data/yam_cdc_hi_egg_titers.tsv', 'data/yam_crick_hi_titers.tsv']
            with open(out, 'w+') as f:
                call = ['cat'] + hi_titers_yam
                print call
                subprocess.call(call, stdout=f)

        # Copy TSVs to nextflu/augur directory. Leave in fauna/data/ for nextstrain/augur.
        if args.copy_data:
            os.system("cp data/h3n2_crick_hi_strains.tsv %sdata/h3n2_hi_strains.tsv"%(args.augur_path))
            os.system("cp data/h3n2_crick_hi_titers.tsv %sdata/h3n2_hi_titers.tsv"%(args.augur_path))
            os.system("cp data/h1n1pdm_crick_hi_strains.tsv %sdata/h1n1pdm_hi_strains.tsv"%(args.augur_path))
            os.system("cp data/h1n1pdm_crick_hi_titers.tsv %sdata/h1n1pdm_hi_titers.tsv"%(args.augur_path))
            os.system("cp data/vic_crick_hi_strains.tsv %sdata/vic_hi_strains.tsv"%(args.augur_path))
            os.system("cp data/vic_crick_hi_titers.tsv %sdata/vic_hi_titers.tsv"%(args.augur_path))
            os.system("cp data/yam_crick_hi_strains.tsv %sdata/yam_hi_strains.tsv"%(args.augur_path))
            os.system("cp data/yam_crick_hi_titers.tsv %sdata/yam_hi_titers.tsv"%(args.augur_path))

    elif args.virus == "ebola":

        call = "python vdb/ebola_download.py -db vdb -v ebola --fstem ebola"
        os.system(call)

        if args.copy_data:
            os.system("cp data/ebola* %sdata/"%(args.augur_path))

    elif args.virus == "dengue":

        call = "python vdb/dengue_download.py -db vdb -v dengue --fstem dengue"
        os.system(call)

        if args.copy_data:
            os.system("cp data/dengue* %sdata/"%(args.augur_path))

    elif args.virus == "zika":

        call = "python vdb/zika_download.py -db vdb -v zika --fstem zika"
        os.system(call)

        if args.copy_data:
            os.system("cp data/zika* %sdata/"%(args.augur_path))

    elif args.virus == "h7n9":

        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:PB2 --fstem h7n9_pb2")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:PB1 --fstem h7n9_pb1")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:PA --fstem h7n9_pa")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:HA --fstem h7n9_ha")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:NP --fstem h7n9_np")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:NA --fstem h7n9_na")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:MP --fstem h7n9_mp")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:NS --fstem h7n9_ns")

        if args.copy_data:
            os.system("cp data/h7n9* %sdata/"%(args.augur_path))

    else:
        print("%s is an invalid virus type.\nValid viruses are flu, ebola, dengue, zika, and h7n9."%(args.virus))
