## Download data from Genbank



Start from the existing accession list from Dudas et al. 2018: https://github.com/nextstrain/mers/blob/master/data/MERS_CoV_274_CDS.dates.txt

Call this file `mers_accessions.txt` and put it in the `fauna/data/` directory.

## Upload to fauna

`python2 vdb/mers_upload.py -db test_vdb -v mers --ftype accession --source genbank --locus genome --fname mers_accessions.txt`

## Download from fauna

`python2 vdb/mers_download.py -db test_vdb -v mers --fstem mers --resolve_method choose_genbank`
