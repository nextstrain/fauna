# DENGUE Pipeline Notes

## Upload documents to VDB
1. Download sequences from [LANL](https://hfv.lanl.gov/components/sequence/HCV/search/searchi.html)
  * Select GB Submission Date `>=MM/YYYY` of last VDB update or "Last GenBank Update" in the upper right corner (whichever is earlier).
  * Set the rest of the parameters as shown:  
![Parameters](./download_instructions.png)  
  * Hit "Search"  
  * Select "Save Background Info" and check the box for "Click here to include the sequence."
  ![Save](./download_instructions2.png)  
2. Move downloaded file to `fauna/data`
3. Upload to vdb database
  * `python vdb/zika_upload.py -db vdb -v zika --source vipr --locus genome --fname results.tbl`