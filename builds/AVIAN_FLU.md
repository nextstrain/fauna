# Avian flu pipeline notes

## VDB

### Upload documents to VDB

1. Download sequences and meta information from [GISAID](http://platform.gisaid.org/)
  * In EPIFLU, select for either H7N9 sequences or H5N1 sequences, select `HA` as required segment, select Submission Date >= last upload date to vdb
  * Download at most 5000 isolates at a time, may have to split downloads by submission date
  * Download Isolates as XLS with YYYY-MM-DD date format
  * Download Isolates as "Sequences (DNA) as FASTA"
    * Select all DNA
    * Fasta Header as 0: DNA Accession no., 1: Isolate name, 2: Isolate ID, 3: Segment, 4: Passage details/history, 5: Submitting lab
    * `DNA Accession no. | Isolate name | Isolate ID | Segment | Passage details/history | Submitting lab`
2. Move files to `fauna/data` as `gisaid_epiflu.xls` and `gisaid_epiflu.fasta`.
3. Upload to vdb database
  * `python2 vdb/avian_flu_upload.py -db vdb -v avian_flu --source gisaid --fname gisaid_epiflu`
  * Recommend running with `--preview` to confirm strain names and locations are correctly parsed before uploading
  	* Can add to [geo_synonyms file](source-data/geo_synonyms.tsv) and [flu_fix_location_label file](source-data/flu_fix_location_label.tsv) to fix some of the formatting.

1. Download sequences from [IRD](https://www.fludb.org)
  * Search for Sequences and strains
  * Select Data Type as Strain
  * Enter either "H5N1" or "H7N9" under Subtype
  * Click Search
  * Click download all
  ...
  * Download as `GenomicFastaResults.fasta`

### Download documents from VDB

```bash
python vdb/avian_flu_download.py -db vdb -v avian_flu --select locus:HA subtype:h7n9 --fstem h7n9_ha
```
