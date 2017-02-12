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

# Copy FASTAs to nextflu directory
os.system("cp data/h3n2.fasta ../nextflu/augur/data/")
os.system("cp data/h1n1pdm.fasta ../nextflu/augur/data/")
os.system("cp data/vic.fasta ../nextflu/augur/data/")
os.system("cp data/yam.fasta ../nextflu/augur/data/")

# Copy FASTAs to nextflu-cdc directory
os.system("cp data/h3n2.fasta ../nextflu-cdc/augur/data/")
os.system("cp data/h1n1pdm.fasta ../nextflu-cdc/augur/data/")
os.system("cp data/vic.fasta ../nextflu-cdc/augur/data/")
os.system("cp data/yam.fasta ../nextflu-cdc/augur/data/")

# Copy FASTAs to nextflu-cdc-fra directory
os.system("cp data/h3n2.fasta ../nextflu-cdc-fra/augur/data/")
os.system("cp data/h1n1pdm.fasta ../nextflu-cdc-fra/augur/data/")
os.system("cp data/vic.fasta ../nextflu-cdc-fra/augur/data/")
os.system("cp data/yam.fasta ../nextflu-cdc-fra/augur/data/")

# Download public titers from database
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype h3n2 --select assay_type:hi")
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype h1n1pdm --select assay_type:hi")
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype vic --select assay_type:hi")
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype yam --select assay_type:hi")

# Copy titer tables to nextflu directory
os.system("cp data/h3n2_hi_strains.tsv ../nextflu/augur/data/")
os.system("cp data/h3n2_hi_titers.tsv ../nextflu/augur/data/")
os.system("cp data/h1n1pdm_hi_strains.tsv ../nextflu/augur/data/")
os.system("cp data/h1n1pdm_hi_titers.tsv ../nextflu/augur/data/")
os.system("cp data/vic_hi_strains.tsv ../nextflu/augur/data/")
os.system("cp data/vic_hi_titers.tsv ../nextflu/augur/data/")
os.system("cp data/yam_hi_strains.tsv ../nextflu/augur/data/")
os.system("cp data/yam_hi_titers.tsv ../nextflu/augur/data/")

# Download CDC HI titers from database
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype h3n2 --select assay_type:hi")
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype h1n1pdm --select assay_type:hi")
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype vic --select assay_type:hi")
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype yam --select assay_type:hi")

# Copy titers to nextflu-cdc directory
os.system("cp data/h3n2_hi_strains.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/h3n2_hi_titers.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/h1n1pdm_hi_strains.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/h1n1pdm_hi_titers.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/vic_hi_strains.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/vic_hi_titers.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/yam_hi_strains.tsv ../nextflu-cdc/augur/data/")
os.system("cp data/yam_hi_titers.tsv ../nextflu-cdc/augur/data/")

# Download CDC FRA titers from database
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype h3n2 --select assay_type:fra")
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype h1n1pdm --select assay_type:fra")
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype vic --select assay_type:fra")
os.system("python tdb/download.py -db cdc_tdb -v flu --ftype augur --subtype yam --select assay_type:fra")

# Copy titers to nextflu-cdc-fra directory
os.system("cp data/h3n2_hi_strains.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/h3n2_hi_titers.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/h1n1pdm_hi_strains.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/h1n1pdm_hi_titers.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/vic_hi_strains.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/vic_hi_titers.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/yam_hi_strains.tsv ../nextflu-cdc-fra/augur/data/")
os.system("cp data/yam_hi_titers.tsv ../nextflu-cdc-fra/augur/data/")
