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
parser.add_argument('--virus', default="flu", help="virus to download; default is flu")
parser.add_argument('--flu_lineages', default=["h3n2", "h1n1pdm", "vic", "yam"], nargs='+', type = str,  help ="seasonal flu lineages to download, options are h3n2, h1n1pdm, vic and yam")
parser.add_argument('--sequences', default=False, action="store_true", help="download sequences from vdb")
parser.add_argument('--titers', default=False, action="store_true", help="download titers from tdb")
parser.add_argument('--titers_sources', default=["crick", "cdc"], nargs='+', type = str,  help ="titer sources to download, options are crick and cdc")
parser.add_argument('--titers_passages', default=["egg", "cell"], nargs='+', type = str,  help ="titer passage types to download, options are egg and cell")


if __name__=="__main__":
    params = parser.parse_args()

    if params.virus == "flu":
        # Download FASTAs from database
        if params.sequences:
            for lineage in params.flu_lineages:
                call = "python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_%s --fstem %s"%(lineage, lineage)
                print(call)
                os.system(call)

        if params.titers:
            # download titers
            for source in params.titers_sources:
                if source == "crick":
                    for lineage in params.flu_lineages:
                        call = "python tdb/download.py -db tdb -v flu --subtype %s --select assay_type:hi --fstem %s_crick_hi_cell"%(lineage, lineage)
                        print(call)
                        os.system(call)
                if source == "cdc":
                    for passage in params.titers_passages:
                        for lineage in params.flu_lineages:
                            call = "python tdb/download.py -db cdc_tdb -v flu --subtype %s --select assay_type:hi serum_passage_category:%s --fstem %s_cdc_hi_%s"%(lineage, passage, lineage, passage)
                            print(call)
                            os.system(call)
                        lineage = 'h3n2'
                        call = "python tdb/download.py -db cdc_tdb -v flu --subtype %s --select assay_type:fra serum_passage_category:%s --fstem %s_cdc_fra_%s"%(lineage, passage, lineage, passage)
                        print(call)
                        os.system(call)

            # concatenate to create default HI strain TSVs for each subtype
            # need to sum counts across HI strain TSVs
            for lineage in params.flu_lineages:
                strain_to_count = {}
                for source in params.titers_sources:
                    hi_strain_file = 'data/%s_%s_hi_cell_strains.tsv'%(lineage, source)
                    with open(hi_strain_file) as f:
                        lines = f.read().splitlines()
                        for line in lines:
                            (strain, count) = line.split("\t")
                            if strain in strain_to_count:
                                strain_to_count[strain] = strain_to_count[strain] + count
                            else:
                                strain_to_count[strain] = count
                out_file = 'data/%s_hi_strains.tsv'%(lineage)
                with open(out_file, 'w') as f:
                    [f.write('{0}\t{1}\n'.format(strain, count)) for strain, count in strain_to_count.items()]

            # concatenate to create default HI titer TSVs for each subtype
            for lineage in params.flu_lineages:
                out = 'data/%s_hi_titers.tsv'%(lineage)
                hi_titers = []
                for source in params.titers_sources:
                    hi_titers.append('data/%s_%s_hi_cell_titers.tsv'%(lineage, source))
                with open(out, 'w+') as f:
                    call = ['cat'] + hi_titers
                    print call
                    subprocess.call(call, stdout=f)

    elif params.virus == "ebola":

        call = "python vdb/ebola_download.py -db vdb -v ebola --fstem ebola"
        print(call)
        os.system(call)

    elif params.virus == "dengue":

        # Download all serotypes together.
        call = "python vdb/dengue_download.py"
        print(call)
        os.system(call)

        # Download individual serotypes.
        serotypes = [1, 2, 3, 4]
        for serotype in serotypes:
            call = "python vdb/dengue_download.py --select serotype:%i" % serotype
            print(call)
            os.system(call)

        # Download titers.
        if params.titers:
            call = "python tdb/download.py -db tdb -v dengue --fstem dengue"
            print(call)
            os.system(call)

    elif params.virus == "zika":

        call = "python vdb/zika_download.py -db vdb -v zika --fstem zika"
        print(call)
        os.system(call)

    elif params.virus == "mumps":

        call = "python vdb/mumps_download.py -db vdb -v mumps --fstem mumps --resolve_method choose_genbank"
        print(call)
        os.system(call)

    elif params.virus == "h7n9" or params.virus == "avian":

        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:PB2 --fstem h7n9_pb2")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:PB1 --fstem h7n9_pb1")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:PA --fstem h7n9_pa")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:HA --fstem h7n9_ha")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:NP --fstem h7n9_np")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:NA --fstem h7n9_na")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:MP --fstem h7n9_mp")
        os.system("python vdb/h7n9_download.py -db vdb -v h7n9 --select locus:NS --fstem h7n9_ns")

    else:
        print("%s is an invalid virus type.\nValid viruses are flu, ebola, dengue, zika, mumps, h7n9, and avian."%(params.virus))
        sys.exit(2)
