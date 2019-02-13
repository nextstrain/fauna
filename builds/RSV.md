# RSV Pipeline Notes

## Update

* Update citation fields
  * `python2 vdb/update.py -db vdb -v rsv --update_citations`
  * updates `authors`, `title`, `url`, `journal` and `puburl` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again

## Download

    python2 vdb/download.py -db vdb -v rsv --fstem rsv --resolve_method choose_genbank

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=SubmitForm&blockId=1337&decorator=paramyxo)

1. Select all genomes
2. Download sequences
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: Genbank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python2 vdb/rsv_upload.py -db vdb -v rsv --source genbank --locus genome --fname GenomeFastaResults.fasta`
