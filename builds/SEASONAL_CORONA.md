## Download data from Vipr

* [Vipr search URL](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ModifySearch&selectionContext=1584731753653)
* This is search fields of `genome[data to return] AND human[host] under Coronaviridae`
* Select all results. Download as `Genome FASTA` and select `Custom Format`, then add all possible FASTA fields in the order: `genbank_accession|strain_name|segment|date|host|country|subtype|virus_species`
* Download the file as `human_cov_genome.fasta` and save it in the fauna/data/ directory
* By default, `seasonal_corona_upload.py` will try to upload full length sequencing data (the `human_cov_genome.fasta` file downloaded from Vipr) AND gene-only sequences for Spike and HE (obtained through alignment of full-length sequences). Download these files from the [seasonal-corona repo](https://www.github.com/nextstrain/seasonal-corona/starter_data) and save them in the fauna/data/ directory
* If you DO NOT wish to upload gene-only sequences, add `--vipr_upload` to the command below. If you wish to ONLY upload gene sequences (not full length) add `--gene_upload` to the command below 

## Upload to fauna

`python2 vdb/seasonal_corona_upload.py -db vdb -v seasonal_corona --source genbank --fname human_cov_annotated.fasta`

## Download from fauna

`python2 vdb/seasonal_corona_download.py -db vdb -v seasonal_corona --fstem seasonal_corona`
