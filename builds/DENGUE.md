# DENGUE Pipeline Notes

## Upload documents to VDB
1. Download sequences from [LANL](https://hfv.lanl.gov/components/sequence/HCV/search/searchi.html)
  * Select GB Submission Date `>=MM/YYYY` of last VDB update or "Last GenBank Update" in the upper right corner (whichever is earlier).
  * Set the rest of the parameters as shown:  
![Parameters](figures/download_instructions.png)  
  * Hit "Search"  
  * Select "Save Background Info" and check the box for "Click here to include the sequence."  
2. Move downloaded file to `fauna/data`
3. Upload to vdb database
  * `python2 vdb/dengue_upload.py -db vdb -v dengue --fname results.tbl --ftype tsv`

## Download sequence documents from VDB

* `python2 vdb/dengue_download.py` # all serotypes together
* `python2 vdb/dengue_download.py --select serotype:1` # just serotype 1
* `python2 vdb/dengue_download.py --select serotype:2` # just serotype 2
* `python2 vdb/dengue_download.py --select serotype:3` # just serotype 3
* `python2 vdb/dengue_download.py --select serotype:4` # just serotype 4

## Download titer documents from TDB

* `python2 tdb/download.py -db tdb -v dengue --fstem dengue`
