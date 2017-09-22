# Mumps Pipeline Notes

## Obtaining / preparing sequences

#### Publically available sequences (ViPR)
* https://www.viprbrc.org/brc/home.spg?decorator=vipr
* Search -> sequences & strains -> Paramyxoviridae
* settings: genome, search for "mumps" and wait until "Mumps virus (Species)" pops up (103 genomes as of sept 2017), min genome length 14kb
* select all -> download
* download genome fasta with custom format. Fields: GB-SN-SEG-DATE-HOST-COUN-SUBTY-VIRSPEC
* copy `GenomeFastaResults.fasta` to `fauna/data/mumps_vipr.fasta`
* upload to fauna
  * `python vdb/mumps_upload.py -db vdb -v mumps --source genbank --locus genome --fname mumps_vipr.fasta --fasta_header_fix source-data/mumps_header_fix.tsv`

#### Private data in CSV format
