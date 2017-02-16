# Simple script to run required operations to
# 1. Download FASTAs from database
# 2. Copy FASTAs to nextflu directory
# 3. Download titer tables from database
# 4. Copy titer tables to nextflu directory
# Run from base fauna directory with python flu/download_all.py
# Assumes that nextflu/, nextflu-cdc/ and nextflu-cdc-fra/ are
# sister directories to fauna/

import os

# Download FASTAs from database
os.system("python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_h3n2 --fstem h3n2")
os.system("python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_h1n1pdm --fstem h1n1pdm")
os.system("python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_vic --fstem vic")
os.system("python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_yam --fstem yam")

# Copy FASTAs to nextflu/augur directory. Leave in fauna/data/ for nextstrain/augur.
os.system("cp data/h3n2.fasta ../nextflu/augur/data/")
os.system("cp data/h1n1pdm.fasta ../nextflu/augur/data/")
os.system("cp data/vic.fasta ../nextflu/augur/data/")
os.system("cp data/yam.fasta ../nextflu/augur/data/")

# Download Crick titers from database
os.system("python tdb/download.py -db tdb -v flu --subtype h3n2 --select assay_type:hi --fstem h3n2_crick_hi")
os.system("python tdb/download.py -db tdb -v flu --subtype h1n1pdm --select assay_type:hi --fstem h1n1pdm_crick_hi")
os.system("python tdb/download.py -db tdb -v flu --subtype vic --select assay_type:hi --fstem vic_crick_hi")
os.system("python tdb/download.py -db tdb -v flu --subtype yam --select assay_type:hi --fstem yam_crick_hi")

# Copy TSVs to nextflu/augur directory. Leave in fauna/data/ for nextstrain/augur.
os.system("cp data/h3n2_crick_hi_strains.tsv ../nextflu/augur/data/h3n2_hi_strains.tsv")
os.system("cp data/h3n2_crick_hi_titers.tsv ../nextflu/augur/data/h3n2_hi_titers.tsv")
os.system("cp data/h1n1pdm_crick_hi_strains.tsv ../nextflu/augur/data/h1n1pdm_hi_strains.tsv")
os.system("cp data/h1n1pdm_crick_hi_titers.tsv ../nextflu/augur/data/h1n1pdm_hi_titers.tsv")
os.system("cp data/vic_crick_hi_strains.tsv ../nextflu/augur/data/vic_hi_strains.tsv")
os.system("cp data/vic_crick_hi_titers.tsv ../nextflu/augur/data/vic_hi_titers.tsv")
os.system("cp data/yam_crick_hi_strains.tsv ../nextflu/augur/data/yam_hi_strains.tsv")
os.system("cp data/yam_crick_hi_titers.tsv ../nextflu/augur/data/yam_hi_titers.tsv")

# Download CDC HI titers from database
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype h3n2 --select assay_type:hi --fstem h3n2_cdc_hi")
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype h1n1pdm --select assay_type:hi --fstem h1n1pdm_cdc_hi")
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype vic --select assay_type:hi --fstem vic_cdc_hi")
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype yam --select assay_type:hi --fstem yam_cdc_hi")

# Download CDC FRA titers from database
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype h3n2 --select assay_type:fra --fstem h3n2_cdc_fra")
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype h1n1pdm --select assay_type:fra --fstem h1n1pdm_cdc_fra")
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype vic --select assay_type:fra --fstem vic_cdc_fra")
os.system("python tdb/download.py -db cdc_tdb -v flu --subtype yam --select assay_type:fra --fstem yam_cdc_fra")
