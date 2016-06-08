# Notes for ZIBRA project

## Schema notes

* `strain`: This study / strain ID in the form of `ZBN47`, `ZBJP1`, etc... This is the *primary key* of the table and is required for every document.
* `amplicon_concentration`: Purity of DNA after PCR amplification. Measured in ng/ul.
* `citations`: Not crucial. This can be populated with citation information.
* `country`: All samples should be `brazil`.
* `ct`: Ct value of positive RT-PCR result.
* `date`: Collection date of the sample. Should be formatted as `2015-07-27` (YYYY-MM-DD). If a sample lacks complete date information, enter as `2015-07-XX` (day unknown) or `2015-XX-XX` (month and day unknown).
* `division`: Name of Brazilian state of patient origin. Should be snakecase without accents. Examples:
 * `bahia`
 * `para`
 * `paraiba`
 * `rio_grande_do_norte`
* `host`: Host species. Human samples are `human`.
* `lacen_id`: This is LACEN sample ID. Should be a 12-digit number. This very important to be able to link sample to clinical metadata.
* `location`: Name of municipality of patient origin. Should be snakecase without accents. Examples:
 * `canguaretama`
 * `natal` 
 * `nova_cruz`
 * `sao_goncalo_do_amarante`
* `microcephaly`: Whether sample was linked to microcephaly, `true` or `false`.
* `minion_barcode`: List of MinION library/barcodes associated with sample, for example `[2_NB08, 2_NB09]`.
* `onset_date`: Date of symptom onset. Should be formatted as `2015-07-27` (YYYY-MM-DD). 
* `patient_age`: Patient age, `30y`, `6m`, `10d`, etc...
* `patient_sex`: Patient sex, `male` or `female`.
* `public`: Whether the sample genome can be shared publicly, `true` or `false`.
* `region`: All samples should be `south_america`.
* `rt_positive`: Whether RT-PCR is positive for Zika, `true` or `false`.
* `sequences`: Contains three fields:
 * `accession`: Same as `strain`, based on LACEN sample ID.
 * `locus`: All samples should be `genome`.
 * `sequence`: Genome sequence.
* `timestamp`: The last edit time of the database entry in `2016-06-04-14-06` (YYYY-MM-DD-HH-MM) format. This should be
automatically managed by scripts/chateau.
* `virus`: All samples should be `zika`.

Modifications to the schema should be made to [zibra_metadata_upload.py]() and [zibra_download.py]().

## Sample IDs

* `ZBRA1` is sample 1 from location A (LACEN Natal).
* `ZBRB1` is sample 1 from location B (LACEN Joao Pessoa).

## Database commands

Remember to [install rethinkdb bindings](README.md#install).

Upload metadata with:

    python vdb/zibra_metadata_upload.py -db vdb -tb zibra --fname lacen_rn.tsv --ftype tsv --source zibra --virus zika --country brazil --authors ZiBRA --local

Upload sequences with:

    python vdb/zibra_upload.py -db vdb -tb zibra --fname minion.fasta --ftype fasta --source zibra --virus zika --locus genome --local

Download metadata with:

    python vdb/zibra_download.py -db vdb -tb zibra --fstem zibra --ftype tsv --local
    
Download just metadata for samples from `natal`:

    python vdb/zibra_download.py -db vdb -tb zibra --fstem zibra --ftype tsv --select location:natal --local

Push local rethinkdb `vdb.zibra` documents to remote `vdb.zibra` rethinkdb table:
	
	python vdb/sync.py --push --local_table vdb.zibra --remote_table vdb.zibra
	
Pull remote rethinkdb `vdb.zibra` documents to local `vdb.zibra` rethinkdb table:

	python vdb/sync.py --pull --local_table vdb.zibra --remote_table vdb.zibra

## Download latest metadata for consensus builds

Remember to [install rethinkdb bindings](README.md#install).

From `nextstrain-db/` run:

    source environment_rethink.sh
    python vdb/zibra_download.py -db vdb -tb zibra --fstem zibra --ftype tsv

This will result in the file `nextstrain-db/data/zibra.tsv` that has all necessary metadata. This file can be searched for `2_NB07`, etc... in the `minion_barcode` column to match MinION output to metadata, including strain name.
