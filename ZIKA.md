# ZIKA Pipeline Notes

## Upload documents to VDB
1. Download sequences from [VIPR](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_zika)
  * Select year >= 2013 and genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `nextstrain-db/data`
3. Upload to vdb database
  * `python vdb/zika_upload.py -db vdb -v zika --source vipr --locus genome --fname GenomeFastaResults.fasta`
  

## Update documents in VDB
* Update citation fields
  * `python vdb/zika_update.py -db vdb -v zika --update_citations`
  * updates `authors`, `title` and `url` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again
* Update location fields
  * After hand editing `location` in [chateau](https://github.com/blab/chateau)
  * `python vdb/zika_update.py -db vdb -v zika --update_locations`
  * Updates `division`, `country`, `region`, `latitude`, `longitude` fields
  
## Download documents from VDB
* `python vdb/zika_download.py -db vdb -v zika --fstem zika`
