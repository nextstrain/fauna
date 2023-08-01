# Avian flu pipeline notes

## VDB

### Upload documents to VDB

#### Upload from GISAID

1. Download sequences and meta information from [GISAID](http://platform.gisaid.org/)
  * In EPIFLU, select for H9N2, H7N9, or H5 sequences, select `HA` as required segment, select Submission Date >= last upload date to vdb
  * Download at most 5000 isolates at a time, may have to split downloads by submission date
  * Download Isolates as XLS with YYYY-MM-DD date format
  * Download Isolates as "Sequences (DNA) as FASTA"
    * Select all DNA, except HE and PE
    * Fasta Header as 0: DNA Accession no., 1: Isolate name, 2: Isolate ID, 3: Segment, 4: Passage details/history, 5: Submitting lab
    * `DNA Accession no. | Isolate name | Isolate ID | Segment | Passage details/history  | DNA INSDC`
2. Move files to `fauna/data` as `gisaid_epiflu.xls` and `gisaid_epiflu.fasta`.
3. Upload to vdb database
  * `python3 vdb/avian_flu_upload.py -db vdb -v avian_flu --data_source gisaid --source gisaid --fname gisaid_epiflu`
  * Recommend running with `--preview` to confirm strain names and locations are correctly parsed before uploading
  	* Can add to [geo_synonyms file](source-data/geo_synonyms.tsv) and [flu_fix_location_label file](source-data/flu_fix_location_label.tsv) to fix some of the formatting.

#### Upload from IRD

1. Download sequences from [IRD](https://www.fludb.org)
  * Search for Sequences and strains
  * Select Data Type as Strain
  * Enter either "H5N1" or "H7N9" under Subtype
  * Click Search
  * Click download all
  ...
  * Download "Segment FASTA" as `GenomicFastaResults.fasta`. Select "Custom format", select all and add.
2. Move file to `fauna/data` as `GenomicFastaResults.fasta`.
3. Upload to vdb database
  * `python3 vdb/avian_flu_upload.py -db vdb -v avian_flu --data_source ird --source ird --fname GenomicFastaResults.fasta`
  * Recommend running with `--preview` to confirm strain names and locations are correctly parsed before uploading
  	* Can add to [geo_synonyms file](source-data/geo_synonyms.tsv) and [flu_fix_location_label file](source-data/flu_fix_location_label.tsv) to fix some of the formatting.

### Download documents from VDB

```
python3 vdb/avian_flu_download.py -db vdb -v avian_flu --select locus:HA subtype:h7n9 --fstem h7n9_ha
```
