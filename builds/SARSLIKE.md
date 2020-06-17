# ZIKA Pipeline Notes

## Download

    python2 vdb/download.py -db vdb -v coronavirus --fstem coronavirus --resolve_method choose_genbank

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=corona)

#### For all beta-CoV

1. Download sequences
  * Select genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python2 vdb/coronavirus_upload.py -db vdb -v coronavirus --source genbank --locus genome --fname GenomicFastaResults.fasta`

#### For SARS-like CoV

1. Download sequences
  * Select genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python2 vdb/coronavirus_upload.py -db vdb -v sarslike --source genbank --locus genome --fname GenomicFastaResults.fasta`  

### [Wuhan genome from Virological](http://virological.org/t/initial-genome-release-of-novel-coronavirus/319)

Upload with:

    python2 vdb/coronavirus_upload.py -db vdb -v sarslike --source virological --locus genome --authors "Zhang et al" --url http://virological.org/t/initial-genome-release-of-novel-coronavirus/319 --title "Initial genome release of novel coronavirus" --fname WH-Human_1.fasta

### [5 Wuhan genomes from GISAID](https://gisaid.org)

Upload with:

    python2 vdb/coronavirus_upload.py -db vdb -v sarslike --source gisaid --locus genome --url https://www.gisaid.org/ --title "Newly discovered betacoronavirus, Wuhan 2019-2020" --fname sarslike_gisaid.fasta

## Update

* Update citation fields
  * `python2 vdb/coronavirus_update.py -db vdb -v coronavirus --update_citations`
