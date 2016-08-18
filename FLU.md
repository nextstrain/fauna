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
2. Rename files to same name with `.xls` and `.fasta` and move to `nextstrain-db/vdb/data`
  * ie. `gisaid_2002-2004.xls` and `gisaid_2002-2004.fasta`
  * `--fname` for upload will be the base name ie. `gisaid_2002-2004`
3. Upload to vdb database
  * `python vdb/flu_upload.py -db vdb -v flu --source gisaid --fname gisaid_2002-2004`
  * Recommend running with `--preview` to confirm strain names and locations are correctly parsed before uploading
  	* Can add to [geo_synonyms file](source-data/geo_synonyms.tsv), [flu_strain_name_fix file](source-data/flu_strain_name_fix.tsv) and [flu_fix_location_label file](source-data/flu_fix_location_label.tsv) to fix some of the formatting.

### Update documents in VDB
* Update genetic grouping fields
  * `python vdb/flu_update.py -db vdb -v flu --update_groupings`
  * updates `vtype`, `subtype`, `lineage`
  
### Download documents from VDB
* `python vdb/zika_download.py -db vdb -v zika`

## TDB

### Upload documents to TDB
1. Convert [NIMR report](https://www.crick.ac.uk/research/worldwide-influenza-centre/annual-and-interim-reports/) pdfs to csv files
2. Move csv files to subtype directory in `nextstrain-db/tdb/data/`
3. Upload to tdb database
  * `python tdb/upload.py -db tdb -v flu --subtype h1n1pdm`
  * Recommend running with `--preview` to confirm strain names are correctly parsed before uploading
  	* Can add to [HI_ref_name_abbreviations file](source-data/HI_ref_name_abbreviations.tsv) and [HI_flu_strain_name_fix file](source-data/HI_flu_strain_name_fix.tsv) to fix some strain names.

### Download documents from TDB
* `python tdb/download.py -db tdb -v flu --subtype h1n1pdm`
* Downloads all titer measurements and counts of HI_strains
