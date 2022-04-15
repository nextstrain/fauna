# ZIKA Pipeline Notes

## Setup

1. Make sure environment variables for connecting to fauna are set.

## Upload via ViPR and update citations

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_zika)

1. Download sequences
  * Select year >= 2013 and genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
  * May also use the [ViPR API](https://www.viprbrc.org/brc/staticContent.spg?decorator=reo&type=ViprInfo&subtype=API)

  ```
  curl "https://www.viprbrc.org/brc/api/sequence?datatype=genome&family=flavi&species=Zika%20virus&fromyear=2013&minlength=5000&metadata=genbank,strainname,segment,date,host,country,genotype,species&output=fasta" |\
  tr '-' '_' |\
  tr ' ' '_' |\
  sed 's:N/A:NA:g' >\
  GenomicFastaResults.fasta
  ```
  
  The search-and-replace commands (`tr`, `sed`) are necessary because the API downloads fasta headers similar to:

  `>KY241742|ZIKV_SG_072|N/A|2016-08-28|Human|Singapore|Asian|Zika virus`
  
  but need to match the GUI downloaded headers similar to: 
  
  `>KY241742|ZIKV_SG_072|NA|2016_08_28|Human|Singapore|Asian|Zika_virus`


2. Move downloaded sequences to `fauna/data`
3. Extract `GenomicFastaResults.tar.gz` and rename the extracted file to `GenomicFastaResults.fasta`
4. Upload to vdb database
  * `python3 vdb/zika_upload.py -db vdb -v zika --source genbank --locus genome --fname GenomicFastaResults.fasta`

### Update

* Update citation fields
  * `python3 vdb/zika_update.py -db vdb -v zika --update_citations`
  * updates `authors`, `title`, `url`, `journal` and `puburl` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again

## Download from Fauna, parse, compress and push to S3

### Download from Fauna

```
python3 vdb/download.py \
  --database vdb \
  --virus zika \
  --fasta_fields strain virus accession collection_date region country division location source locus authors url title journal puburl \
  --resolve_method choose_genbank \
  --fstem zika
```

This results in the file `data/zika.fasta` with FASTA header ordered as above.

### Parse

```
augur parse \
  --sequences data/zika.fasta \
  --output-sequences data/sequences.fasta \
  --output-metadata data/metadata.tsv \
  --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  --prettify-fields region country division city
```

This results in the files `data/sequences.fasta` and `data/metadata.tsv`.

### Compress

```
xz --compress data/sequences.fasta
gzip data/metadata.tsv
```

This results in the files `data/sequences.fasta.xz` and `data/metadata.tsv.gz`.

### Push to S3

```
nextstrain remote upload s3://nextstrain-data/files/zika/ data/sequences.fasta.xz data/metadata.tsv.gz
```

This pushes files to S3 to be made available at https://data.nextstrain.org/files/zika/sequences.fasta.xz and https://data.nextstrain.org/files/zika/metadata.tsv.gz.

## Run zika workflow

See instructions at https://github.com/nextstrain/zika.
