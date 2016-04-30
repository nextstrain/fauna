# TDB
The titer database (TDB) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database. 

## Uploading
Titer measurements from the Crick Worldwide Influenza Centre can be uploaded from csv files that are in a standard format.

## Schema

Each document in the database represents an HI test between a virus and serum from a specific ferret. Multiple titer measurements are stored if the test is repeated.

* `virus`: strain name of virus tested.
* `serum`: strain name serum was raised against. 
* `ferret_id`: id of the ferret the serum was raised in.
* `source`: List of document names from which titer measurements were added. 
* `passage`: Passage history of the virus.
* `index`: Used as compound index in rethinkdb. Array with `virus`, `serum`, `ferret_id`, `source`, `passage`.
* `titer`: List of all titer measurements for this test. 
* `date_modified`:  Last modification date for document in `YYYY-MM-DD` format.
* `date`: Collection date of virus in `YYYY-MM-DD` format, for example, `2016-02-28`.
* `group`: Genetic group of the virus. 
* `ref`: Boolean for whether the virus is a reference virus. `True` if reference virus. 

### Attribute Requirements
Tests with null values for required attributes will be filtered out of those uploaded. Viruses with missing optional attributes will still be uploaded
* Required attributes: `virus`, `serum`, `ferret_id`, `index`, `titer`, `source`, `date_modified`
* Optional attributes: `date`, `passage`, `group`, `ref`

### Commands
Command line arguments to run vdb_upload:
* -db --database default='tdb', help=database to upload to. Ex 'tdb', 'test'
* -v --virus help=virus table to interact with. Ex 'h1n1pdm', 'h3n2', 'vic', 'yam'
* --overwrite default=False help=whether to overwrite existing non-null fields
* --path help=path to fasta file, default is data/virus/
* --auth\_key help=authorization key for rethink database
* --host help=rethink host url

### Examples:

python tdb/src/tdb_upload.py -db tdb -v h1n1pdm

## Downloading
Measurements can be downloaded from tdb
* Downloads all documents in database
* Prints result to designated fasta or json file. 
	* Writes null attributes as '?'
	* Writes text file description in this order (0:`virus`, 1:`serum`, 2:`ferret_id`, 3:`source`, 4:`titer`)

### Commands
Command line arguments to run vdb_download:
* -db --database default='vdb', help=database to download from. Ex 'vdb', 'test'
* -v --virus help=virus table to interact with. Ex 'Zika', 'Flu'
* --path help=path to dump output files to, default is data/
* --ftype help=output file format, default is 'fasta', other option is 'json'
* --fstem help=output file stem name, default is VirusName\_Year\_Month\_Date
* --auth\_key help=authorization key for rethink database
* --host help=rethink host url

### Examples:

python tdb/src/tdb_download.py -db tdb -v H1N1pdm

## Accessing the Database
All titer measurements are stored using [Rethinkdb deployed on AWS](https://www.rethinkdb.com/docs/paas/#deploying-on-aws)

To access tdb you need an authorization key. This can be passed as a command line argument (see above) or set as an environment variable with a bash script.

`source environment_rethink.sh`
```shell
#!/bin/bash
export RETHINK_AUTH_KEY=EXAMPLE_KEY
export RETHINK_HOST=EXAMPLE_HOST
export NCBI_EMAIL=example\@email.org
```