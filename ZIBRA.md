# Notes for ZIBRA project

## Excel "schema"

One sample per row. Missing or inapplicable fields should be left blank. Highlighting is okay, but do not include any textual notes except for in the `notes` field. Suggest blue for positive RT and orange for microcephaly.

* `strain`: This study / strain ID in the form of `ZBRA1`, `ZBRB1`, etc... This is the *primary key* of the table and is required for every sample.
* `sample_number`: This is the sample number from this collection location, `1`, `2`, `3`, etc...
* `lab`: Collection lab, `lacen_natal`, `lacen_joao_pessoa`, `fiocruz_recife`, etc...
* `lab_id`: This is LACEN sample ID (or other institution sample ID). Should be a 12-digit number. This very important to be able to link sample to clinical metadata. Each LACEN uses its own IDs.
* `rt_positive`: Whether RT-PCR is positive for Zika, `true` or `false`.
* `ct`: Ct value of positive RT-PCR result.
* `extraction_date`: Date of RNA extraction.
* `amplicon_concentration`: Purity of DNA after PCR amplification. Measured in ng/ul.
* `minion_barcodes`: List of MinION library/barcodes associated with sample separated by commas, for example `2_nb08,2_nb09`.
* `country`: All samples should be `brazil`.
* `state`: Name of Brazilian state of patient origin. Should be snakecase without accents, `bahia`, `rio_grande_do_norte`, etc...
* `municipality`: Name of municipality of patient origin. Should be snakecase without accents, `canguaretama`, `natal`, `sao_goncalo_do_amarante`, etc...
* `collection_date`: Collection date of the sample. Should be formatted as `2015-07-27` (YYYY-MM-DD). If a sample lacks complete date information, enter as `2015-07-XX` (day unknown) or `2015-XX-XX` (month and day unknown).
* `onset_date`: Date of symptom onset. Should be formatted as `2015-07-27` (YYYY-MM-DD).
* `host_species`: Host species. Human samples are `human`.
* `patient_age`: Patient age, `30y`, `6m`, `10d`, etc...
* `patient_sex`: Patient sex, `male` or `female`.
* `pregnant`: Whether the patient was pregnant, `true` or `false`.
* `pregnancy_week`: Week of pregnancy, as number `10`, `12`, etc....
* `pregnancy_trimester`: Trimester of pregnancy `1`, `2`, `3`.
* `microcephaly`: Whether sample was linked to microcephaly, `true` or `false`.
* `sample_type`: Type of sample, `serum`, `urine`, etc...
* `symptoms`: Free-form text of clinical symptoms.
* `notes`: Free-form text for notes associated with this sample.

Modifications to the schema should be made to [zibra_metadata_upload.py]() and [zibra_download.py]().

## Sample IDs

* `ZBRA1` is sample 1 from location A (LACEN Natal).
* `ZBRB1` is sample 1 from location B (LACEN Joao Pessoa).

## Database commands

Remember to [install rethinkdb bindings](README.md#install).

Upload metadata with:

    python vdb/zibra_metadata_upload.py -db vdb -tb zibra --fname zibra.tsv --ftype tsv --source zibra --virus zika --country brazil --authors ZiBRA --local

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
