# VDB
The virus database (VDB) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

## Uploading

Sequences can be uploaded from a fasta file, genbank file or file of genbank accession number to a virus specific table within vdb. It currently:

* Uploads from an input file
	* FASTA
		* Reads fasta description in this order for `zika_upload` (0:`accession`, 2:`strain`, 4:`date`, 6:`country`)
		* Reading order can easily be changed in source code for different viruses or files
	* Genbank
		* Can include multiple genbank entries
		* Parses information from entries
	* Accession
		* One accession number on each line in file
		* Uses entrez to get genbank entries from accession numbers
* Uploads from command line argument `--accessions`
	* comma separated list of accession numbers
	* Uses entrez to get genbank entries from accession numbers
* Formats information to fit with vdb schema.
* Uploads information to virus table
	* Appends to list of sequences if new accession number. If no accession number, appends if new sequence.
	* If strain already in database, default is to update attributes with new information only if current attribute is null. Can overwrite existing non-null data with `--overwrite` option.

## Schema

Information about the virus is stored in one table (ie. flu\_viruses) with information about corresponding
sequence information stored in another table (ie. flu\_sequences). They are linked together by the `sequences`
field in the viruses table that stores a list of accession numbers of corresponding sequences.

Flu Virus Table:

* `strain`: primary key. The canonical strain name. For flu this would be something like `A/Perth/16/2009`.
* `isolate_id`: used for some viruses, like those from gisaid. `EPI_ISL_221430`.
* `virus`: virus type in CamelCase format. Loose term for like viruses (viruses that you'd want to include in a single tree). Examples include `flu`, `ebola`, `zika`.
* `subtype`: virus subtype in lowercase, where available, Null otherwise. `h3n2`, `h1n1pdm`, `vic`, `yam`
* `collection_date`: collection date in `YYYY-MM-DD` format, for example, `2016-02-28` or `2016-02-XX` if day ambiguous.
* `region`: collection region where virus was isolated in snakecase format.  See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_regions.tsv) for examples.
* `country`: collection country where virus was isolated in snakecase format. See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_synonyms.tsv) for examples.
* `division`: administrative division where virus was isolated in snakecase format. Where available, Null otherwise.
* `location`: specific location where virus was isolated in snakecase format. Where available, Null otherwise.
* `originating_lab`: lab where the virus was originally isolated from `alaska_state_virology_lab`
* `gender`: gender of the individual from which the virus was isolated `male`, `female`
* `age`: age of the individual from which the virus was isolated, with the age unit. `47y`
* `host`: species of the individual from which the virus was isolated. `human`, `swine`
* `sequences`: list of sequence accession numbers associated with this virus from the sequences table
* `number_of_sequences`: number of sequences for this virus
* `timestamp`: last modification date and time for virus document in `YYYY-MM-DD` format.

Flu Sequences Table:

* `accession`: primary key. The accession number. Assigned a unique random accession number if the sequence doesn't have one.
* `strain`: the canonical strain name to link back to the viruses table. For flu this would be something like `A/Perth/16/2009`.
* `source`: database source for the sequence`genbank`, `gisaid`, etc...
* `passage`: passage history
* `locus`: gene or genomic region, `ha`, `na`, `genome`, etc...
* `sequence`: actual sequence. Lower case.
* `authors`: authors to attribute credit to `Azevedo et al`
* `title`: title of reference. `Discovery of a persistent Zika virus lineage in Bahia, Brazil`
* `url`: url of reference if available, search crossref database for DOI, or link to genbank entry. `http://dx.doi.org/10.1101/049916`, `http://www.ncbi.nlm.nih.gov/nuccore/ku365779`
* `public`: false if sequence not be be released publicly. True otherwise.
* `submitting_lab`: lab that produced and submitted the sequence `who_national_influenza_centre_russian_federation`
* `submission_date`: submission date in `YYYY-MM-DD` format, for example, `2016-02-28` or `2016-02-XX` if day ambiguous.
* `timestamp`: last modification date and time for virus document in `YYYY-MM-DD` format.

### Commands

Command line arguments to run `upload.py`:

* `-db --database`: database to upload to, eg. `vdb`, `test_vdb`
* `--host`: rethink host url
* `--auth\_key`: authorization key for rethink database
* `--local`: include when connecting to a local rethinkdb instance
* `-v --virus`: name of virus and table eg. `zika`, `flu`
* `--fname`: input file name
* `--ftype`: input file type, fasta, genbank or accession
* `--path`: path to fasta file, default is `data/`
* `--accessions`: comma separated list of accessions numbers to upload
* `--email`: to upload viruses via accession must include email to use entrez
* `--overwrite`: overwrite existing non-null fields
* `--preview`: include to preview documents without uploading
* `--replace`: include to replace all documents in database,

Assign attribute to all viruses and sequences being uploaded with these arguments
* `--source`: database source for the sequence`genbank`, `gisaid`, etc...
* `--locus`: gene or genomic region, `ha`, `na`, `genome`, etc...
* `--authors`: authors of source of sequences
* `--country`: country where virus was isolated
* `--private`: to designate sequences being uploaded as not `public`


### Examples

Upload flu sequences from GISAID:

    python vdb/gisaid_flu_upload.py -db test_vdb -v flu --fname 2016_gisaid --source gisaid

Upload Zika sequences from VIPR:

    python vdb/zika_upload.py -db vdb -v zika --source genbank --locus genome --fname GenomeFastaResults.fasta

Upload via accession file:

	python vdb/zika\_upload.py -db test\_vdb -v zika --ftype accession --source genbank --locus genome --virus zika --fname entrez_test.txt

Upload via accession list:

	python vdb/zika\_upload.py -db test\_vdb --v zika --source genbank --locus genome  --virus zika --accessions KU501216,KU501217,KU365780,KU365777

## Downloading

Sequences can be downloaded from vdb.

* Downloads all documents in database
* Each sequence has associated virus meta data paired with it
* Prints result to designated fasta or json file.
	* Writes null attributes as '?'
	* Writes fasta description in this order by default(0:`strain`, 1:`virus`, 2:`accession`, 3:`date`, 4:`region`, 5:`country`, 6:`division`, 7:`location`, 8:`source`, 9:`locus`, 10:`authors`, 11:`subtype`)

### Commands

Command line arguments to run `download.py`:

* `-db --database`: database to download from, eg. `vdb`, `test_vdb`
* `--host`: rethink host url
* `--auth\_key`: authorization key for rethink database
* `--local`: include when connecting to a local rethinkdb instance
* `-v --virus`, virus table to interact with, eg. `zika`, `flu`
* `--ftype`: output file format, default is `fasta`, other option is `json`
* `--fstem`: output file stem name, default is `VirusName\_Year\_Month\_Date`
* `--path`: path to dump output files to, default is `data/`

Subset documents with these commands
* `--public\_only`: include to subset public sequences
* `--select`: Select specific fields to be certain values eg. `--select field1:value1 field2:value1,value2`
* `--present`: Select specific fields to be non-null eg. `--present field1 field2`
* `--interval`: Select date fields to be in a certain interval eg. `--interval collection_date:2016-01-01,2016-01-15`

Prevent the default of resolving duplicate sequences for the same locus with thi argument
* `--keep_duplicates`: keep all duplicate sequences for the same locus

### Examples

Download sequences for `Zika_process.py`:

    python vdb/zika_download.py -db vdb -v zika --fstem zika

    python vdb/zika_download.py -db vdb -v zika --fstem zika --countries brazil haiti --public_only

Download sequences from `flu_download.py`:

    python vdb/download.py -db vdb -v zika --fstem zika

    python vdb/download.py -db vdb -v zika --ftype json --countries brazil haiti --public_only

    python vdb/flu_download.py -db test_vdb -v flu --select locus:HA --present age --interval collection_date:2016-01-01,2016-01-15

## Updating

Sequences and Viruses in vdb can be updated.
* Default update
  * Finds all sequences in current database whose source is `genbank` or `vipr`
  * Uses [entrez](http://www.ncbi.nlm.nih.gov/books/NBK25501/) to update virus and sequence documents
  * `python vdb/zika_update.py -db vdb -v zika`
* `--update_citations`
  * updates `authors`, `title`, `url` fields from genbank files
  * If you get `ERROR: Couldn't connect with entrez, please run again` just run command again
  * `python vdb/zika_update.py -db vdb -v zika --update_citations`
* `--update_locations`
  * First manually edit most detailed location field (ie `location`) with [chateau](https://github.com/blab/chateau)
  * Updates `division`, `country`, `region`, `latitude`, `longitude` fields
  * `python vdb/zika_update.py -db vdb -v zika --update_locations`
* `--update_groupings`
  * Updates genetic grouping fields like `vtype`, `subtype`, `lineage
  * Only implemented currently for [flu_update.py](vdb/flu_update.py)
  * `python vdb/flu_update.py -db vdb -v flu --update_groupings`

### Examples

	python vdb/zika_update.py -db vdb -v zika

	python vdb/update.py -db test_vdb -v zika --accessions KU501216,KU501217

## Backup and Restore

VDB tables can be backed up to S3 or locally.

* Backups can be run manually or continuously everyday
* Backs up all tables in database
* Restoration keeps current documents in database but overwrites conflicting documents with the same primary key


### Examples

Backup every `vdb` table to s3 backup file

	python vdb/backup.py -db vdb --backup_s3

Restore `vdb.zika` to 2016-05-25 version from s3 backup file

	python vdb/restore.py -db vdb -v zika --backup_s3 --restore_date 2016-05-25

Backup every `vdb` table to local backup file

	python vdb/backup.py -db vdb --backup_local

Restore `vdb.zika` to 2016-05-25 version from local backup file

	python vdb/restore.py -db vdb -v zika --backup_local --restore_date 2016-05-25

Backup `vdb` to s3 backup file everyday

	python vdb/backup.py -db vdb --continuous_backup --backup_s3

## Append

VDB documents can be appended to other tables.

### Examples

Append `vdb` zika documents to `test_vdb` zika tables:

	python vdb/append.py -v zika --from_database vdb --to_database test_vdb


## Sync

VDB tables can be synced between a local rethinkdb instance and external rethinkdb instance.

Push local rethinkdb `test_vdb.zika` documents to remote `vdb.zika` rethinkdb table:

	python vdb/sync.py --push --local_table test_vdb.zika --remote_table test_vdb.zika

Pull remote rethinkdb `test_vdb.zika` documents to local `test_vdb.zika` rethinkdb table:

	python vdb/sync.py --pull --local_table test_vdb.zika --remote_table test_vdb.zika

## Accessing the Database

All viruses are stored using [Rethinkdb deployed on AWS](https://www.rethinkdb.com/docs/paas/#deploying-on-aws). To access vdb you need an authorization key. This can be passed as a command line argument (see above) or set as an environment variable with a bash script, by running `source environment_rethink.sh`:

```shell
#!/bin/bash
export RETHINK_AUTH_KEY=EXAMPLE_KEY
export RETHINK_HOST=EXAMPLE_HOST
export NCBI_EMAIL=example\@email.org
```
