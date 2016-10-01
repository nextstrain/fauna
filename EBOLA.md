# EBOLA Pipeline Notes

## Upload documents to VDB

1. Using current set of Ebola sequences from ebola.nextstrain.org.
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python vdb/ebola_upload.py -db vdb -v ebola --source various --locus genome --fname ebola_transfer.fasta`

## Download documents from VDB
  * `python vdb/ebola_download.py -db vdb -v ebola --fstem ebola`
