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
    * `DNA Accession no. | Isolate name | Isolate ID | Segment | Passage details/history | Submitting lab`
2. Move files to `fauna/data` as `gisaid_epiflu.xls` and `gisaid_epiflu.fasta`.
3. Upload to vdb database
  * `python vdb/flu_upload.py -db vdb -v flu --source gisaid --fname gisaid_epiflu`
  * Recommend running with `--preview` to confirm strain names and locations are correctly parsed before uploading
  	* Can add to [geo_synonyms file](source-data/geo_synonyms.tsv), [flu_strain_name_fix file](source-data/flu_strain_name_fix.tsv) and [flu_fix_location_label file](source-data/flu_fix_location_label.tsv) to fix some of the formatting.

### Update documents in VDB

All of these functions are quite slow given they run over ~600k documents. Use sparingly.

* Update genetic grouping fields
  * `python vdb/flu_update.py -db vdb -v flu --update_groupings`
  * updates `vtype`, `subtype`, `lineage`

* Update locations
  * `python vdb/flu_update.py -db vdb -v flu --update_locations`
  * updates `division`, `country` and `region` from `location`

* Update passage_category fields
  * `python vdb/flu_update.py -db vdb -v flu --update_passage_categories`
  * update `passage_category` based on `passage` field

### Download documents from VDB

* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_h3n2 --fstem h3n2`
* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_h1n1pdm --fstem h1n1pdm`
* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_vic --fstem vic`
* `python vdb/flu_download.py -db vdb -v flu --select locus:HA lineage:seasonal_yam --fstem yam`

## TDB

### Upload documents to TDB

#### Raw tables from NIMR reports

1. Convert [NIMR report](https://www.crick.ac.uk/research/worldwide-influenza-centre/annual-and-interim-reports/) pdfs to csv files
2. Move csv files to subtype directory in `fauna/data/`
3. Upload to tdb database
  * `python tdb/upload.py -db tdb -v flu --subtype h3n2 --ftype flat --fstem h3n2_nimr_titers`
  * Recommend running with `--preview` to confirm strain names are correctly parsed before uploading
  	* Can add to [HI_ref_name_abbreviations file](source-data/HI_ref_name_abbreviations.tsv) and [HI_flu_strain_name_fix file](source-data/HI_flu_strain_name_fix.tsv) to fix some strain names.

#### Flat files

1. Move line-list tsv files to `fauna/data/`
2. Upload to tdb database with `python tdb/upload.py -db tdb -v flu --subtype h3n2 --ftype flat --fstem H3N2_HI_titers_upload`

#### CDC files

1. Move line-list tsv files to `fauna/data/`
2. Upload HI titers to tdb database with `python tdb/cdc_upload.py -db cdc_tdb -v flu --ftype flat --fstem HITest_Oct2016_to_Sep2017_titers`
3. Upload FRA titers to tdb database with `python tdb/cdc_upload.py -db cdc_tdb -v flu --ftype flat --fstem FRA_Oct2016_to_Sep2017_titers`

#### Crick files

1. Move Excel documents to `fauna/data/`
2. Run `python tdb/crick_upload.py -db crick_tdb --assay_type hi --fstem H3N2HIs`
3. Run `python tdb/crick_upload.py -db crick_tdb --assay_type fra --fstem H3N2VNs`
4. Run `python tdb/crick_upload.py -db crick_tdb --assay_type hi --fstem H1N1pdm09HIs`
5. Run `python tdb/crick_upload.py -db crick_tdb --assay_type hi --fstem BVicHIs`
6. Run `python tdb/crick_upload.py -db crick_tdb --assay_type hi --fstem BYamHIs`

#### NIID files

1. Make sure `NIID-Tokyo-WHO-CC/` is a sister directory to `fauna/`
2. Upload all titers with `python tdb/upload_all.py --sources niid -db niid_tdb`

#### VIDRL files

1. Make sure `VIDRL-Melbourne-WHO-CC/` is a sister directory to `fauna/`
2. Upload all titers with `python tdb/upload_all.py --sources vidrl -db vidrl_tdb`

### Download documents from TDB

* `python tdb/download.py -db tdb -v flu --subtype h3n2`
* `python tdb/download.py -db tdb -v flu --subtype h1n1pdm`
* `python tdb/download.py -db tdb -v flu --subtype vic`
* `python tdb/download.py -db tdb -v flu --subtype yam`
