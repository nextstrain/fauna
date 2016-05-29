# VDB
The virus database (VDB) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

## Uploading
Sequences can be uploaded from a fasta file, genbank file or file of genbank accession number to a virus specific table within vdb. It currently
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

* `Strain`: primary key. The canonical strain name. For flu this would be something like `A/Perth/16/2009`.
* `Virus`: Virus type in CamelCase format. Loose term for like viruses (viruses that you'd want to include in a single tree). Examples include `flu`, `ebola`, `zika`.
* `Subtype`: Virus subtype in lowercase, where available, Null otherwise. `h3n2`, `h1n1pdm`, `vic`, `yam`
* `Date_Modified`: Last modification date for virus document in `YYYY-MM-DD` format.
* `Date`: Collection date in `YYYY-MM-DD` format, for example, `2016-02-28` or `2016-02-xx` if day ambiguous.
* `Region`: Collection region in CamelCase format.  See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_regions.tsv) for examples. 
* `Country`: Collection country in CamelCase format. See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_synonyms.tsv) for examples.
* `Division`: Administrative division in CamelCase format. Where available, Null otherwise.
* `Location`: Specific location in CamelCase format. Where available, Null otherwise.
* `Public`: True if publicly available sequence. False otherwise.
* `Sequences`: list of sequences...
  * `Accession`: Accession number. Where available, Null otherwise.
  * `Source`: Genbank, GISAID, etc... in CamelCase format.
  * `Locus`: gene or genomic region, `HA`, `NA`, `Genome`, etc... in CamelCase format.
  * `Sequence`: Actual sequence. Upper case.
*`Citations`: list of citations for corresponding sequences
  * `Authors`: Authors to attribute credit to. Where available, Null otherwise. in CamelCase format.
  * `Title`: Title of reference.
  * `url`: Url of reference if available, search crossref database for DOI, otherwise link to genbank entry. 

### Attribute Requirements
Viruses with null values for required attributes will be filtered out of those uploaded. Viruses with missing optional attributes will still be uploaded
* Required virus attributes: `strain`, `date`, `country`, `sequences`, `virus`, `date_modified`, `public`
* Required sequence attributes: `source`, `locus`, `sequence`
* Optional virus attributes: `division`, `location`
* Optional sequence attributes: `accession`, `authors`, `title`, `url`

### Commands
Command line arguments to run vdb_upload:
* -db --database default='vdb', help=database to upload to. Ex 'vdb', 'test'
* -v --virus help=virus table to interact with. Ex 'zika', 'zlu'
* --fname help=input file name
* --ftype help=input file type, fasta, genbank or accession
* --accessions help=comma separated list of accessions numbers to upload
* --source
* --locus
* --authors help=authors of source of sequences
* --private help=to designate sequences being uploaded as not `public`
* --overwrite default=False help=whether to overwrite existing non-null fields
* --path help=path to fasta file, default is data/virus/
* --auth\_key help=authorization key for rethink database
* --host help=rethink host url
* --email help=to upload viruses via accession must include email to use entrez

### Examples:

Upload flu sequences from GISAID:

    python vdb/flu_upload.py -db test_vdb -v flu --fname gisaid_epiflu_sequence.fasta --source gisaid

Upload Zika sequences from VIPR:

    python vdb/zika_upload.py --database vdb --virus zika --fname GenomeFastaResults.fasta --source genbank --locus genome
    
Upload via accession file:

	python vdb/zika_upload.py --database test --virus zika --fname entrez_test.txt --ftype accession --source genbank --locus genome

Upload via accession list:

	python vdb/zika_upload.py --database test --virus zika --source genbank --locus genome --accessions KU501216,KU501217,KU365780,KU365777

## Downloading
Sequences can be downloaded from vdb.
* Downloads all documents in database
* If virus has more than one sequence, picks the longest sequence
* Prints result to designated fasta or json file. 
	* Writes null attributes as '?'
	* Writes fasta description in this order (0:`strain`, 1:`virus`, 2:`accession`, 3:`date`, 4:`region`, 5:`country`, 6:`division`, 7:`location`, 8:`source`, 9:`locus`, 10:`authors`, 11:`subtype`)

###Commands
Command line arguments to run vdb_download:
* -db --database default='vdb', help=database to download from. Ex 'vdb', 'test'
* -v --virus help=virus table to interact with. Ex 'zika', 'flu'
* --path help=path to dump output files to, default is data/
* --ftype help=output file format, default is 'fasta', other option is 'json'
* --fstem help=output file stem name, default is VirusName\_Year\_Month\_Date
* --auth\_key help=authorization key for rethink database
* --host help=rethink host url
* --public\_only help=include to subset public sequences
* --countries help=Countries(in CamelCase Format) to be include in download, multiple arguments allowed

### Examples:

Download sequences for `Zika_process.py`:

    python vdb/download.py -db vdb -v zika --fstem zika
    
    python vdb/download.py -db test -v zika --ftype json --countries Brazil Haiti --public_only

## Updating
Sequences in vdb can be automatically updated
* Only sequences whose source is Genbank
* Uses entrez to check for updates to certain fields
* Updates the fields: `authors`, `title`, `url`, `sequence` 
* Must specify database, virus and email.

### Examples:

	python vdb/update.py -db test -v zika
	
	python vdb/update.py -db test -v zika --accessions KU501216,KU501217
	
## Backup and Restore
VDB tables can be backed up to S3 or to a local source
* Can be run manually or continuously everyday
* Deletes old backups after a certain length (default is 50 days)
* Restoration keeps current documents in database

### Examples

Backup `test_vdb.zika` to s3 backup file
	
	python vdb/backup.py -db test_vdb --backup --backup_s3

Restore `test_vdb.zika` to 2016-05-25 from s3 backup file
	
	python vdb/backup.py -db test_vdb --restore --restore_table zika --restore_date 2016-05-25 --backup_s3

Backup `test_vdb.zika` to local backup file
	
	python vdb/backup.py -db test_vdb --backup --backup_local

Restore `test_vdb.zika` to 2016-05-25 from local backup file
	
	python vdb/backup.py -db test_vdb --restore --restore_table zika --restore_date 2016-05-25 --backup_local

Backup `test_vdb.zika` to s3 backup file everyday	
	
	python vdb/backup.py -db test_vdb --continuous_backup --backup_s3
	
## Append
VDB documents can be appended to other tables

### Examples
	
Append `vdb.zika` documents to `test_vdb.zika`

	python vdb/append.py --from_table vdb.zika --to_table test_vdb.zika

## Sync
VDB tables can be synced between a local rethinkdb instance and external rethinkdb instance

Push local rethinkdb test_vdb.zika documents to remote vdb.zika rethinkdb table
	
	python vdb/sync.py --push --local_table test_vdb.zika --remote_table test_vdb.zika	

Pull remote rethinkdb test_vdb.zika documents to local test_vdb.zika rethinkdb table

	python vdb/sync.py --pull --local_table test_vdb.zika --remote_table test_vdb.zika

## Accessing the Database
All viruses are stored using [Rethinkdb deployed on AWS](https://www.rethinkdb.com/docs/paas/#deploying-on-aws)

To access vdb you need an authorization key. This can be passed as a command line argument (see above) or set as an environment variable with a bash script.

`source environment_rethink.sh`
```shell
#!/bin/bash
export RETHINK_AUTH_KEY=EXAMPLE_KEY
export RETHINK_HOST=EXAMPLE_HOST
export NCBI_EMAIL=example\@email.org
```

