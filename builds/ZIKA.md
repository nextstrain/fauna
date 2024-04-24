# ZIKA Pipeline Notes

_The use of Zika in Fauna is now defunct. Instead of uploading Zika sequences and metadata to Fauna, it's directly downloaded from NCBI via `nextstrain/zika/ingest` and uploaded to S3 as part of the ingest pipeline. The notes below are for future reference, but are not actively maintained or acted upon._

## Ingest data from NCBI GenBank

Navigate to the nextstrain/zika repository and [follow the instructions for ingest](https://github.com/nextstrain/zika/tree/main/ingest).

```
git clone https://github.com/nextstrain/zika.git
cd zika
cd ingest
nextstrain build .
```

This results in the files `results/metadata.tsv` and `results/sequences.fasta`

## Compress

```
zstd -T0 results/sequences.fasta
zstd -T0 results/metadata.tsv
```

This results in the files `results/sequences.fasta.zst` and `results/metadata.tsv.zst`.

## Upload data to s3

```
nextstrain remote upload s3://nextstrain-data/files/workflows/zika/ results/sequences.fasta.zst
nextstrain remote upload s3://nextstrain-data/files/workflows/zika/ results/metadata.tsv.zst
```

This pushes files to S3 to be made available at https://data.nextstrain.org/files/workflows/zika/sequences.fasta.zst and https://data.nextstrain.org/files/workflows/zika/metadata.tsv.zst.

## Run zika workflow

See instructions at https://github.com/nextstrain/zika/tree/main/phylogenetic

```
cd ../phylogenetic
mv ingest/results data
nextstrain build .
```
