# ZIKA Pipeline Notes

## Download

    python2 vdb/download.py -db vdb -v ncov --fstem ncov --resolve_method choose_genbank

## Upload

### [Genomes from GISAID](https://gisaid.org)

Upload with:

    python2 vdb/ncov_upload.py -db vdb -v ncov --source gisaid --locus genome --url https://www.gisaid.org --title "Newly discovered betacoronavirus, BetaCoV 2019-2020" --fname ncov_gisaid.fasta

### Single genomes from Genbank

    python2 vdb/ncov_upload.py -db vdb -v ncov --source genbank --locus genome --ftype accession --fname ncov_accessions.txt
