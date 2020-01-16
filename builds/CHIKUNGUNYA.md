# ZIKA Pipeline Notes

## Download

    python2 vdb/download.py -db vdb -v chikungunya --fstem chikungunya --resolve_method choose_genbank

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=SubmitForm&blockId=728&decorator=toga#)

1. Download sequences
  * Select all genomes
  * Download as Genome FASTA
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python2 vdb/chikungunya_upload.py -db vdb -v chikungunya --source genbank --locus genome --fname GenomicFastaResults.fasta`

## Update

* Update citation fields
  * `python2 vdb/zika_update.py -db vdb -v zika --update_citations`
  * updates `authors`, `title`, `url`, `journal` and `puburl` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again
