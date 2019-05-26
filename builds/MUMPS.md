## Download data from Genbank

* [Genbank search URL](https://www.ncbi.nlm.nih.gov/nuccore?term=mumps%5Btitle%5D%20AND%20viruses%5Bfilter%5D%20AND%20%28%225000%22%5BSLEN%5D%20%3A%20%2220000%22%5BSLEN%5D%29&cmd=DetailsSearch)
* This is search fields of `mumps[title] AND viruses[filter] AND ("5000"[SLEN] : "20000"[SLEN])`
* Send to : Complete Record : File : Accession List
* This downloads the file `sequence.seq`
* Open this file and remove the `.1`, `.2`, etc... from the accession numbers

## Upload to fauna

`python2 vdb/mumps_upload.py -db vdb -v mumps --ftype accession --source genbank --locus genome --fname sequence.seq`


FASTA header field ordering:
1. random numbering - this will later be filled in by GenBank accession
2. strain name
3. collection date
4. host species
5. country
6. state/region
7. genotype

## Update fauna database

_This is not necessary when uploading accessions as we do here._
This is needed to populate certain attributes such as author & paper title.
`python2 vdb/mumps_update.py -db vdb -v mumps --update_citations`

## Download from fauna

`python2 vdb/mumps_download.py -db vdb -v mumps --fstem mumps --resolve_method choose_genbank`

## Upload Broad genomes

Preprocess to fix metadata and header ordering

`python2 vdb/mumps_preprocess_fasta.py --fasta data/muv-nextstrain-20170718.pruned.fasta > data/mumps_broad.fasta`

Upload to fauna

`python2 vdb/mumps_upload.py -db vdb -v mumps --source broad --locus genome --fname mumps_broad.fasta --authors "Wohl et al" --title "Unpublished"`

## Upload BCCDC genomes

If you have a FASTA file and CSV metadata, this script will help (with minor modifications as needed)

`python2 scripts/mumps.csv-and-fasta-to-vipr-fasta.py data/input.mumps.raw.fasta data/input.mumps.csv data/input.mumps.vipr.fasta`


Upload to fauna

`python2 vdb/mumps_upload.py -db vdb -v mumps --source bccdc --locus genome --fname mumps.bc.fasta --authors "Gardy et al" --title "Unpublished"`

## Upload Fred Hutch genomes

Upload to fauna

`python2 vdb/mumps_upload.py -db vdb -v mumps --source fh --locus genome --fname MuVs-WA0268502_buccal-Washington.USA-16.fasta --authors "Moncla et al" --title "Unpublished"`
