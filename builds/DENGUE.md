# DENGUE Pipeline Notes

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_dengue)

1. Download sequences
  * Select genome length >= 5000
  * Download as Genome FASTA
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Type
2. Move downloaded sequences to `fauna/GenomicFastaResults.fasta`
3. Upload to vdb database
  * `python3 vdb/dengue_upload.py -db vdb -v dengue --source genbank --locus genome --fname GenomicFastaResults.fasta`

## Update

* Update citation fields
  * `python3 vdb/dengue_update.py -db vdb -v dengue --update_citations`
  * updates `authors`, `title`, `url`, `journal` and `puburl` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again    

## Download sequence documents from VDB

* `python3 vdb/dengue_download.py -v dengue` # all serotypes together
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_1` # just serotype 1
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_2` # just serotype 2
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_3` # just serotype 3
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_4` # just serotype 4

## Download titer documents from TDB

* `python3 tdb/download.py -db tdb -v dengue --fstem dengue`
