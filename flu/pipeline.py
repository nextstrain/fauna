# simple script to run required operations to
# 1. Upload GISAID XLS and FASTA to database
# 2. Update genetic groupings
# 3. Download FASTAs from database
# 4. Copy FASTAs to nextflu directory
# 5. Download titer tables from database
# 6. Copy titer tables to nextflu directory
# Run from base nextstrain-db directory with python flu/pipeline.py
# Requires files gisaid_epiflu.xls and gisaid_epiflu.fasta in data/

import os

# Upload GISAID XLS and FASTA to database
os.system("python vdb/flu_upload.py -db vdb -v flu --source gisaid --fname gisaid_epiflu")

# Update genetic groupings
os.system("python vdb/flu_update.py -db vdb -v flu --update_groupings")

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

# Download titer tables from database
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype h3n2")
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype h1n1pdm")
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype vic")
os.system("python tdb/download.py -db tdb -v flu --ftype augur --subtype yam")

# Copy titer tables to nextflu directory
os.system("cp data/h3n2_hi_strains.txt ../nextflu/augur/data/")
os.system("cp data/h3n2_hi_titers.txt ../nextflu/augur/data/")
os.system("cp data/h1n1pdm_hi_strains.txt ../nextflu/augur/data/")
os.system("cp data/h1n1pdm_hi_titers.txt ../nextflu/augur/data/")
os.system("cp data/vic_hi_strains.txt ../nextflu/augur/data/")
os.system("cp data/vic_hi_titers.txt ../nextflu/augur/data/")
os.system("cp data/yam_hi_strains.txt ../nextflu/augur/data/")
os.system("cp data/yam_hi_titers.txt ../nextflu/augur/data/")
