# Simple script to run required operations to
# 1. Download FASTAs from database
# 2. Copy FASTAs to nextflu directory
# 3. Download titer tables from database
# 4. Copy titer tables to nextflu directory
# Run from base fauna directory with python flu/download_all.py
# Assumes that nextflu/, nextflu-cdc/ and nextflu-cdc-fra/ are
# sister directories to fauna/

import os
copy_data = False
sequences=False
crick=False
cdc=True

# Download FASTAs from database
if sequences:
    for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
        call = "python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_%s --fstem %s"%(lineage, lineage)
        os.system(call)

# Copy FASTAs to nextflu/augur directory. Leave in fauna/data/ for nextstrain/augur.
if copy_data:
    os.system("cp data/h3n2.fasta ../nextflu/augur/data/")
    os.system("cp data/h1n1pdm.fasta ../nextflu/augur/data/")
    os.system("cp data/vic.fasta ../nextflu/augur/data/")
    os.system("cp data/yam.fasta ../nextflu/augur/data/")

if crick:
    for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
        call = "python tdb/download.py -db tdb -v flu --subtype %s --select assay_type:hi --fstem %s_crick_hi"%(lineage, lineage)
        os.system(call)


# Download CDC HI titers from database
if cdc:
    for lineage in ['h3n2', 'h1n1pdm', 'vic', 'yam']:
        for passage in ["egg", "cell"]:
            call = "python tdb/download.py -db cdc_tdb -v flu --subtype %s --select assay_type:hi serum_passage_category:%s --fstem %s_cdc_hi_%s"%(lineage, passage, lineage, passage)
            os.system(call)

    # Download CDC FRA titers from database
    lineage = 'h3n2'
    for passage in ["egg", "cell"]:
        call = "python tdb/download.py -db cdc_tdb -v flu --subtype %s --select assay_type:fra serum_passage_category:%s --fstem %s_cdc_fra_%s"%(lineage, passage, lineage, passage)
        os.system(call)

# Copy TSVs to nextflu/augur directory. Leave in fauna/data/ for nextstrain/augur.
if copy_data:
    os.system("cp data/h3n2_crick_hi_strains.tsv ../nextflu/augur/data/h3n2_hi_strains.tsv")
    os.system("cp data/h3n2_crick_hi_titers.tsv ../nextflu/augur/data/h3n2_hi_titers.tsv")
    os.system("cp data/h1n1pdm_crick_hi_strains.tsv ../nextflu/augur/data/h1n1pdm_hi_strains.tsv")
    os.system("cp data/h1n1pdm_crick_hi_titers.tsv ../nextflu/augur/data/h1n1pdm_hi_titers.tsv")
    os.system("cp data/vic_crick_hi_strains.tsv ../nextflu/augur/data/vic_hi_strains.tsv")
    os.system("cp data/vic_crick_hi_titers.tsv ../nextflu/augur/data/vic_hi_titers.tsv")
    os.system("cp data/yam_crick_hi_strains.tsv ../nextflu/augur/data/yam_hi_strains.tsv")
    os.system("cp data/yam_crick_hi_titers.tsv ../nextflu/augur/data/yam_hi_titers.tsv")
