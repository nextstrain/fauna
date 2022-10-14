## Download data from Genbank

* [Genbank search URL](https://www.ncbi.nlm.nih.gov/nuccore?term=measles%5Btitle%5D%20AND%20viruses%5Bfilter%5D%20AND%20%28%225000%22%5BSLEN%5D%20%3A%20%2220000%22%5BSLEN%5D%29&cmd=DetailsSearch)
* This is search fields of `measles[title] AND viruses[filter] AND ("5000"[SLEN] : "20000"[SLEN])`
* Send to : Complete Record : File : Accession List
* This downloads the file `sequence.seq`
* Remove the `.1`, `.2`, etc... from the accession numbers in `sequence.seq`:

    ```
    sed -i '' -e 's/.1$//g' -e 's/.2$//g' sequence.seq
    ```

## Upload to fauna

```
python3 vdb/measles_upload.py \
  -db vdb \
  -v measles \
  --ftype accession \
  --source genbank \
  --locus genome \
  --fname sequence.seq
```

## Download from fauna

```
python3 vdb/download.py \
  --database vdb \
  --virus measles \
  --fasta_fields strain virus accession collection_date region country division location source locus authors url title journal puburl \
  --resolve_method choose_genbank \
  --fstem measles
```

This results in the file `data/measles.fasta` with FASTA header ordered as above.

### Parse

```
augur parse \
  --sequences data/measles.fasta \
  --output-sequences data/sequences.fasta \
  --output-metadata data/metadata.tsv \
  --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  --prettify-fields region country division city
```

This results in the files `data/sequences.fasta` and `data/metadata.tsv`.

### Compress

```
zstd -T0 data/sequences.fasta
zstd -T0 data/metadata.tsv
```

This results in the files `data/sequences.fasta.zst` and `data/metadata.tsv.zst`.

### Push to S3

```
nextstrain remote upload s3://nextstrain-data/files/measles/ data/sequences.fasta.zst data/metadata.tsv.zst
```

This pushes files to S3 to be made available at https://data.nextstrain.org/files/measles/sequences.fasta.zst and https://data.nextstrain.org/files/measles/metadata.tsv.zst.

## Run measles workflow

See instructions at https://github.com/nextstrain/measles.