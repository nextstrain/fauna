## Download data from Genbank

* [Genbank search URL](https://www.ncbi.nlm.nih.gov/nuccore?term=measles%5Btitle%5D%20AND%20viruses%5Bfilter%5D%20AND%20%28%225000%22%5BSLEN%5D%20%3A%20%2220000%22%5BSLEN%5D%29&cmd=DetailsSearch)
* This is search fields of `measles[title] AND viruses[filter] AND ("5000"[SLEN] : "20000"[SLEN])`
* Send to : Complete Record : File : Accession List
* This downloads the file `sequence.seq`
* Open this file and remove the `.1`, `.2`, etc... from the accession numbers

## Upload to fauna

`python2 vdb/measles_upload.py -db vdb -v measles --ftype accession --source genbank --locus genome --fname sequence.seq`

## Download from fauna

`python2 vdb/measles_download.py -db vdb -v measles --fstem measles --resolve_method choose_genbank`
