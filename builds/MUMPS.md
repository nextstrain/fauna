## Download data from ViPR

_This is deprecated. Instead download data from Genbank below and upload directly._
* Goto [ViPR Paramyxo](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=paramyxo)
* Type "mumps" in search
* Select min length of 5000
* Search
* Select all ->  Download
* Select "Genome FASTA" with custom format. Fields:
  1. Genbank Accession
  2. Strain Name
  3. Segment
  4. Date
  5. Host
  6. Country
  7. Subtype
  8. Virus Species
* Copy `GenomeFastaResults.fasta` file to `fauna/data/GenomeFastaResults.fasta`

## Clean FASTA

_This is deprecated. Instead download data from Genbank below and upload directly._
```
cd fauna
python vdb/mumps_preprocess_fasta.py --vipr
```
The cleaned fasta is now available at `data/mumps_vipr.fasta`

## Download data from Genbank

* [Genbank search URL](https://www.ncbi.nlm.nih.gov/nuccore?term=mumps%5Btitle%5D%20AND%20viruses%5Bfilter%5D%20AND%20%28%225000%22%5BSLEN%5D%20%3A%20%2220000%22%5BSLEN%5D%29&cmd=DetailsSearch)
* This is search fields of `mumps[title] AND viruses[filter] AND ("5000"[SLEN] : "20000"[SLEN])`
* Send to : Complete Record : File : Accession List
* This downloads the file `sequence.seq`
* Open this file and remove the `.1`, `.2`, etc... from the accession numbers

## Upload to fauna

`python vdb/mumps_upload.py -db vdb -v mumps --ftype accession --source genbank --locus genome --fname sequence.seq`

## Update fauna database

_This is not necessary when uploading accessions as we do here._
This is needed to populate certain attributes such as author & paper title.
`python vdb/mumps_update.py -db vdb -v mumps --update_citations`

## Download from fauna

`python vdb/mumps_download.py -db vdb -v mumps --fstem mumps --resolve_method choose_genbank`

## Upload Broad genomes

Preprocess to fix metadata and header ordering

`python vdb/mumps_preprocess_fasta.py --fasta data/muv-nextstrain-20170718.pruned.fasta > data/mumps_broad.fasta`

Upload to fauna

`python vdb/mumps_upload.py -db vdb -v mumps --source broad --locus genome --fname mumps_broad.fasta --authors "Wohl et al" --title "Unpublished"`

## Upload BCCDC genomes

Upload to fauna

`python vdb/mumps_upload.py -db vdb -v mumps --source bccdc --locus genome --fname mumps.bc.fasta --authors "Gardy et al" --title "Unpublished"`
