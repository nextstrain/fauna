# DENGUE Pipeline Notes

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_dengue)

1. Download sequences
  * Select genome length >= 5000
  * Download as Genome FASTA
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Type

3. Move downloaded sequences to `fauna/data/GenomicFastaResults.fasta`

  ```
  cd fauna
  mkdir -p data
  mv ~/Downloads/GenomicFastaResults.tar.gz data/.
  tar -xf data/GenomicFastaResults.tar.gz
  mv data/*GenomicFastaResults.fasta data/GenomicFastaResults.fasta
  ```
  
4. Upload to vdb database

  ```
  python3 vdb/dengue_upload.py \
    -db vdb \
    -v dengue \
    --source genbank \
    --locus genome \
    --fname GenomicFastaResults.fasta
  ```

## Update

* Update citation fields

  ```
  python3 vdb/dengue_update.py \
    -db vdb \
    -v dengue \
    --update_citations
  ```
  
  * updates `authors`, `title`, `url`, `journal` and `puburl` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again

## Download sequence documents from VDB

* `python3 vdb/dengue_download.py -v dengue` # all serotypes together
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_1` # just serotype 1
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_2` # just serotype 2
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_3` # just serotype 3
* `python3 vdb/dengue_download.py -v dengue --select serotype:dengue_virus_4` # just serotype 4

## Download titer documents from TDB

* `python3 tdb/download.py -db tdb -v dengue --fstem dengue`

## Prepare for s3 upload

```
# Download dengue_all
python3 vdb/download.py \
  --database vdb \
  --virus dengue \
  --fasta_fields strain virus accession collection_date region country division location source locus authors url title journal puburl \
  --resolve_method choose_genbank \
  --fstem dengue_all

# Convert to sequences and metadata files
augur parse \
  --sequences data/dengue_all.fasta \
  --output-sequences data/sequences_all.fasta \
  --output-metadata data/metadata_all.tsv \
  --fields strain virus accession date region country division city db segment authors url title journal paper_url \
  --prettify-fields region country division city

# Loop through the 4 serotypes
ARR=(1 2 3 4)

for TYPE in "${ARR[@]}"; do
  python3 vdb/download.py \
    --database vdb \
    --virus dengue \
    --fasta_fields strain virus accession collection_date region country division location source locus authors url title journal puburl \
    --select serotype:Dengue_virus_${TYPE} \
    --fstem dengue_denv${TYPE}

  # Convert to sequences and metadata files
  augur parse \
    --sequences data/dengue_denv${TYPE}.fasta \
    --output-sequences data/sequences_denv${TYPE}.fasta \
    --output-metadata data/metadata_denv${TYPE}.tsv \
    --fields strain virus accession date region country division city db segment authors url title journal paper_url \
    --prettify-fields region country division city
done
```

### Compress

```
ARR=(all denv1 denv2 denv3 denv4)

for SEROTYPE in "${ARR[@]}" do
  zstd -T0 data/sequences_${SEROTYPE}.fasta
  zstd -T0 data/metadata_${SEROTYPE}.tsv
done
```

This results in multiple `data/sequences_*.fasta.zst` and `data/metadata_*.tsv.zst` files.

### Push to S3

```
nextstrain remote upload s3://nextstrain-data/files/dengue/ data/sequences_*.zst data/metadata_*.zst
```

This pushes files to S3 to be made available at https://data.nextstrain.org/files/dengue/sequences_*.fasta.zst and https://data.nextstrain.org/files/dengue/metadata_*.tsv.zst.

## Run dengue workflow

See instructions at https://github.com/nextstrain/dengue.
