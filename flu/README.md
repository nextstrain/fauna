# Flu Pipeline Notes

## VDB

### Upload documents to VDB

1. Download sequences and meta information from [GISAID](http://platform.gisaid.org/)
  * In EPIFLU, select host as `human`, select `HA` as required segment, select Submission Date >= last upload date to vdb
  * Ideally download about 5000 isolates at a time, may have to split downloads by submission date
  * Download Isolates as XLS with YYYY-MM-DD date format
  * Download Isolates as "Sequences (DNA) as FASTA"
    * Select all DNA
    * Fasta Header as 0: DNA Accession no., 1: Isolate name, 2: Isolate ID, 3: Segment, 4: Passage details/history, 5: Submitting lab
2. Move files to `nextstrain-db/data` as `gisaid_epiflu.xls` and `gisaid_epiflu.fasta`.
3. Upload to vdb database
  * `python vdb/flu_upload.py -db vdb -v flu --source gisaid --fname gisaid_epiflu`
  * Recommend running with `--preview` to confirm strain names and locations are correctly parsed before uploading
  	* Can add to [geo_synonyms file](source-data/geo_synonyms.tsv), [flu_strain_name_fix file](source-data/flu_strain_name_fix.tsv) and [flu_fix_location_label file](source-data/flu_fix_location_label.tsv) to fix some of the formatting.

### Update documents in VDB

* Update genetic grouping fields
  * `python vdb/flu_update.py -db vdb -v flu --update_groupings`
  * updates `vtype`, `subtype`, `lineage`
  
### Download documents from VDB

* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_h3n2 --fstem h3n2`
* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_h1n1pdm --fstem h1n1pdm`
* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_vic --fstem vic`
* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_yam --fstem yam`

## TDB

### Upload documents to TDB

1. Convert [NIMR report](https://www.crick.ac.uk/research/worldwide-influenza-centre/annual-and-interim-reports/) pdfs to csv files
2. Move csv files to subtype directory in `nextstrain-db/tdb/data/`
3. Upload to tdb database
  * `python tdb/upload.py -db tdb -v flu --subtype h1n1pdm`
  * Recommend running with `--preview` to confirm strain names are correctly parsed before uploading
  	* Can add to [HI_ref_name_abbreviations file](source-data/HI_ref_name_abbreviations.tsv) and [HI_flu_strain_name_fix file](source-data/HI_flu_strain_name_fix.tsv) to fix some strain names.

### Download documents from TDB

* `python tdb/download.py -db tdb -v flu --ftype augur --subtype h3n2`
* `python tdb/download.py -db tdb -v flu --ftype augur --subtype h1n1pdm`
* `python tdb/download.py -db tdb -v flu --ftype augur --subtype vic`
* `python tdb/download.py -db tdb -v flu --ftype augur --subtype yam`
