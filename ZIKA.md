# ZIKA Pipeline Notes

## Upload documents to VDB
1. Download sequences from [VIPR](https://www.viprbrc.org/brc/vipr_genome_search.spg?method=ShowCleanSearch&decorator=flavi_zika)
  * Select year >= 2013 and genome length >= 5000
  * Download as Genome Fasta
  * Set Custom Format Fields to 0: GenBank Accession, 1: Strain Name, 2: Segment, 3: Date, 4: Host, 5: Country, 6: Subtype, 7: Virus Species
2. Move downloaded sequences to `fauna/data`
3. Upload to vdb database
  * `python vdb/zika_upload.py -db vdb -v zika --source vipr --locus genome --fname GenomeFastaResults.fasta`


## Update documents in VDB
* Update citation fields
  * `python vdb/zika_update.py -db vdb -v zika --update_citations`
  * updates `authors`, `title` and `url` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again
* Update location fields
  * After hand editing `location` in [chateau](https://github.com/blab/chateau)
  * `python vdb/zika_update.py -db vdb -v zika --update_locations`
  * Updates `division`, `country`, `region`, `latitude`, `longitude` fields

## Download documents from VDB
* `python vdb/zika_download.py -db vdb -v zika --fstem zika`

# ZIBRA sequences

## Upload documents to VDB
Regex replace: `^>([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|(\S+) .+`
With: `>\1|\1|\5|brazil|\4|\3`
then
Regex replace: `^>([^|]+)\|\d.+`
With: `>\1|\1|XXXX-XX-XX|brazil`
Upload with: `python vdb/zibra_upload.py -db vdb -v zika --source zibra --locus genome --authors "Zika in Brazil Real-time Analysis Consortium" --fname BRA_ZIBRA_Good.fasta`
And: `python vdb/zibra_upload.py -db vdb -v zika --source zibra --locus genome --authors "Zika in Brazil Real-time Analysis Consortium" --fname BRA_ZIBRA_Partial.fasta`

# Andersen sequences

## Upload documents to VDB
Regex replace: `^>([^|]+)\|XX\|`
With: `>\1|\1|`
Upload with: `python vdb/zibra_upload.py -db vdb -v zika --source andersen --locus genome --authors "Grubaugh et al" --fname andersen.fasta`

# Broad sequences

## Upload documents to VDB
Upload with: `python vdb/zika_upload.py -db vdb -v zika --source broad --locus genome --authors "Broad Viral Genomics Group" --fname ZIKV_BROAD_2016-10.fasta`

# USAMRIID sequences

## Upload documents to VDB
Regex replace: `^>([^_]+)_[^_]+_[^_]+_(\S+)`
With: `>\1|\1|\2|human|usa|florida|florida`
Upload with: `python vdb/zibra_upload.py -db vdb -v zika --source usamriid --locus genome --authors "Ladner et al" --fname RIID_ZIKV_FL_10-27-16.fasta`

# FH sequences

## Upload documents to VDB
Upload with:

* `python vdb/zibra_upload.py -db vdb -v zika --source fh --locus genome --authors "Black et al" --fname ZIKA_good.fasta`
* `python vdb/zibra_upload.py -db vdb -v zika --source fh --locus genome --authors "Black et al" --fname ZIKA_partial.fasta`
