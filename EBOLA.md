# EBOLA Pipeline Notes

## Upload documents to VDB

1. Download Ebola sequence data from https://github.com/ebov/space-time/tree/master/Data
2. Move `Makona_1610_genomes_genomes.fasta` sequences to `fauna/data`
3. Replace `SLE` with `sierra_leone`, `LBR` with `liberia` and `GIN` with `guinea`
4. Replace `sierra_leone\|\?` with `sierra_leone|sierra_leone`, `liberia\|\?` with `liberia|liberia` and `guinea\|\?` with `guinea|guinea`
5. Replace `sierra_leone\|\|` with `sierra_leone|sierra_leone|`, `liberia\|\|` with `liberia|liberia|` and `guinea\|\|` with `guinea|guinea|`
6. Upload to vdb database
  * `python vdb/ebola_upload.py -db vdb -v ebola --source genbank --locus genome --fname Makona_1610_genomes_genbank.fasta`
7. Hand edit author and url info into other 153 genomes
8. Upload to vdb database
  * `python vdb/ebola_upload.py -db vdb -v ebola --source genbank --locus genome --fname Makona_1610_genomes_quick.fasta --authors "Quick et al" --url https://github.com/nickloman/ebov/`

## Update

1. Update citation fields
  * `python vdb/ebola_update.py -db vdb -v ebola --update_citations`
  * Updates `authors`, `title` and `url` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again

## Download documents from VDB
  * `python vdb/ebola_download.py -db vdb -v ebola --fstem ebola --resolve_method choose_genbank`
