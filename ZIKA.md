# ZIKA Pipeline Notes

## Update

* Update citation fields
  * `python vdb/zika_update.py -db vdb -v zika --update_citations`
  * updates `authors`, `title` and `url` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again
* Update location fields
  * After hand editing `location` in [chateau](https://github.com/blab/chateau)
  * `python vdb/zika_update.py -db vdb -v zika --update_locations`
  * Updates `division`, `country`, `region`, `latitude`, `longitude` fields

## Download

    python vdb/zika_download.py -db vdb -v zika --fstem zika

## Upload

### [ViPR sequences](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_zika)

1. Download sequences
  * Select year >= 2013 and genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python vdb/zika_upload.py -db vdb -v zika --source vipr --locus genome --fname GenomeFastaResults.fasta`

### [ZiBRA sequences](https://github.com/zibraproject/zibraproject.github.io/tree/master/data/consensus)

* Replace `^>([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|(\S+) .+` with `>\1|\1|\5|brazil|\4|\3`
* Replace `^>([^|]+)\|\d.+` with `>\1|\1|XXXX-XX-XX|brazil`

Upload with:

    python vdb/zibra_upload.py -db vdb -v zika --source zibra --locus genome --authors "Zika in Brazil Real-time Analysis Consortium" --fname BRA_ZIBRA_Good.fasta
    python vdb/zibra_upload.py -db vdb -v zika --source zibra --locus genome --authors "Zika in Brazil Real-time Analysis Consortium" --fname BRA_ZIBRA_Partial.fasta

### [Scripps sequences](https://github.com/andersen-lab/zika-florida/tree/master/consensus_sequences)

* Replace: `^>([^|]+)\|XX\|` with `>\1|\1|`

Upload with:

    python vdb/zibra_upload.py -db vdb -v zika --source andersen --locus genome --authors "Grubaugh et al" --fname andersen.fasta

### [Broad sequences](http://virological.org/t/33-zika-virus-genomes-sequenced-from-patient-and-pooled-mosquito-samples/372)

Upload with:

    python vdb/zika_upload.py -db vdb -v zika --source broad --locus genome --authors "Broad Viral Genomics Group" --fname ZIKV_BROAD_2016-10.fasta

### [USAMRIID sequences](https://github.com/jtladner/ZIKA_Florida/tree/master/sequences)

* Replace: `^>([^_]+)_[^_]+_[^_]+_(\S+)` with `>\1|\1|\2|human|usa|florida|florida`

Upload with:

    python vdb/zibra_upload.py -db vdb -v zika --source usamriid --locus genome --authors "Ladner et al" --fname RIID_ZIKV_FL_10-27-16.fasta

### [Fred Hutch sequences](https://github.com/blab/zika-seq/tree/master/consensus-genomes)

Upload with:

    python vdb/zibra_upload.py -db vdb -v zika --source fh --locus genome --authors "Black et al" --fname ZIKA_good.fasta
    python vdb/zibra_upload.py -db vdb -v zika --source fh --locus genome --authors "Black et al" --fname ZIKA_partial.fasta
