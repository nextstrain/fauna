# ZIKA Pipeline Notes

## Upload via ViPR and update citations

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_zika)

1. Download sequences
  * Select year >= 2013 and genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python3 vdb/zika_upload.py -db vdb -v zika --source genbank --locus genome --fname GenomicFastaResults.fasta`

### Update

* Update citation fields
  * `python2 vdb/zika_update.py -db vdb -v zika --update_citations`
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

### Compress sequences

```
xz --compress data/sequences.fasta
```

This results in the files `data/sequences.fasta.xz`.

### Push to S3

```
mv data/metadata.tsv data/metadata.tsv.gz
nextstrain remote upload s3://nextstrain-data/files/zika/ data/sequences.fasta.xz data/metadata.tsv.gz
```

This pushes files to S3 to be made available at https://data.nextstrain.org/files/zika/sequences.fasta.xz and https://data.nextstrain.org/files/zika/metadata.tsv.gz. Text files are automatically gzipped on upload with `nextstrain remote upload` and so the file is renamed locally to compensate.
