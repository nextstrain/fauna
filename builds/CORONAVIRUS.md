# ZIKA Pipeline Notes

## Download

    python2 vdb/download.py -db vdb -v coronavirus --fstem coronavirus --resolve_method choose_genbank

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=corona)

1. Download sequences
  * Select genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python2 vdb/coronavirus_upload.py -db vdb -v coronavirus --source genbank --locus genome --fname GenomicFastaResults.fasta`
