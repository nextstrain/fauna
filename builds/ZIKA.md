# ZIKA Pipeline Notes

## Ingest data from NCBI GenBank

Navigate to the nextstrain/zika repository and [follow the instructions for ingest](https://github.com/nextstrain/zika/tree/persephone/ingest).

```
git clone https://github.com/nextstrain/zika.git
cd zika
git checkout persephone
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

> [!NOTE]
> Make sure [authentication for the S3 remote](https://docs.nextstrain.org/projects/cli/en/stable/remotes/s3/#authentication) is configured.

```
nextstrain remote upload s3://nextstrain-data/files/zika/ results/sequences.fasta.zst
nextstrain remote upload s3://nextstrain-data/files/zika/ results/metadata.tsv.zst
```

This pushes files to S3 to be made available at https://data.nextstrain.org/files/zika/sequences.fasta.zst and https://data.nextstrain.org/files/zika/metadata.tsv.zst.

## Run zika workflow

See instructions at https://github.com/nextstrain/zika/tree/persephone/phylogenetic

```
cd ../phylogenetic
mv ingest/results data
nextstrain build .
```
