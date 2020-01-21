# ZIKA Pipeline Notes

## Download

    python2 vdb/download.py -db vdb -v coronavirus --fstem coronavirus --resolve_method choose_genbank

## Upload

### [Genomes from GISAID](https://gisaid.org)

Upload with:

    python2 vdb/ncov_upload.py -db vdb -v ncov --source gisaid --locus genome --url https://www.gisaid.org --title "Newly discovered betacoronavirus, Wuhan 2019-2020" --fname wuhan_gisaid.fasta
