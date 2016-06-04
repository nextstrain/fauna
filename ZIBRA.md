# Notes for ZIBRA project

## Schema notes

* `strain`: This is LACEN sample ID. Should be a 12-digit number.
* `citations`: Not crucial. This can be populated with citation information.
* `country`: All samples should be `brazil`.
* `date`: Collection date of the sample.
* `division`: Name of Brazilian state of patient origin. Should be snakecase. Examples:
 * `bahia`
 * `para`
 * `paraiba`
 * `rio_grande_do_norte`
* `location`: Name of municipality of patient origin. Should be snakecase. Examples:
 * `canguaretama`
 * `natal` 
 * `nova_cruz`
 * `sao_goncalo_do_amarante`
* `public`: `true` or `false`. Whether the sample genome can be shared publicly.
* `region`: All samples should be `south_america`.
* `timestamp`: The last edit time of the database entry in `2016-06-04-14-06` (YYYY-MM-DD-HH-MM) format. This should be automatically managed by scripts/chateau.
* `virus`: All samples should be `zika`.

## Database commands

Upload metadata with:

    python vdb/zibra_metadata_upload.py -db vdb -tb zibra -v zika --fname lacen_rn.tsv --ftype tsv --source zibra --country brazil --local

Push local rethinkdb `vdb.zibra` documents to remote `vdb.zibra` rethinkdb table:
	
	python vdb/sync.py --push --local_table vdb.zibra --remote_table vdb.zibra
	
Pull remote rethinkdb `vdb.zibra` documents to local `vdb.zibra` rethinkdb table:

	python vdb/sync.py --pull --local_table vdb.zibra --remote_table vdb.zibra
