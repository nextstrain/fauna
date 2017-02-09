# TDB

The titer database (TDB) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database.

## Uploading

Titer measurements from [WHO Collaborating Centers](http://www.who.int/influenza/gisrs_laboratory/collaborating_centres/en/) can be uploaded from csv files that are in a standard format.

## Schema

Each document in the database represents an HI test between a virus and serum from a specific ferret. Multiple titer measurements are stored if the test is repeated.

* `virus`: strain name of virus tested.
* `serum`: strain name serum was raised against.
* `ferret_id`: id of the ferret the serum was raised in.
* `source`: List of document names from which titer measurements were added.
* `passage`: Passage history of the virus.
* `subtype`: Subtype of virus, ie h1n1pdm, vic, yam, h3n2
* `host`: Host of virus, ie Human, Swine
* `index`: Used as compound index in rethinkdb. Array with `virus`, `serum`, `ferret_id`, `source`, `passage`, `subtype`, `host`.
* `titer`: List of all titer measurements for this test.
* `date_modified`:  Last modification date for document in `YYYY-MM-DD` format.
* `date`: Collection date of virus in `YYYY-MM-DD` format, for example, `2016-02-28`.
* `ref`: Boolean for whether the virus is a reference virus. `True` if reference virus.

### Attribute Requirements

Tests with null values for required attributes will be filtered out of those uploaded. Viruses with missing optional attributes will still be uploaded.

* Required attributes: `virus`, `serum`, `ferret_id`, `source` , `passage`, `index`, `titer`, `date_modified`
* Optional attributes: `date`, `ref`

### Commands

Command line arguments to run tdb_upload:

* `-db --database`: database to upload to, eg. `tdb`, `test_tdb`
* `-v --virus`: virus table to interact with, eg. `h1n1pdm`, `h3n2`, `vic`, `yam`
* `--overwrite`: overwrite existing non-null fields
* `--exclusive`: downloads all documents in database to see if measurements present, include `--exclusive` to get each document on its own (takes longer, but better if others updating database at same time)
* `--preview`: if included, preview a virus document to be uploaded
* `--replace` if included, delete all documents in table
* `--path`: path to input file, default is `data/virus/`
* `--auth\_key`: authorization key for rethink database
* `--host`: rethink host url

### Examples

Test parsing of HI tables without actually uploading to database:

    python tdb/upload.py -db tdb -v h1n1pdm --preview

Upload measurements to database:

    python tdb/upload.py -db tdb -v h1n1pdm

Replace all measurements in table before uploading:

    python tdb/upload.py -db tdb -v h1n1pdm --replace

### HI Table Troubleshooting

The Francis Crick Institute releases biannual [reports](https://www.crick.ac.uk/research/worldwide-influenza-centre/annual-and-interim-reports/) that include antigenic analyses of the different subtypes of seasonal flu. These tables are in pdf format and must be converted to `csv` format using a pdf converter like [tabular](https://github.com/tabulapdf/tabula) or [okular](https://okular.kde.org/). The reports are not consistent with their column labels and formatting, serum strain names are difficult to parse, and the pdf converters are not perfect. This required manual curation/fixing of the csv files and expansion of the [parse_HI_matrix function](https://github.com/blab/nextflu/blob/master/augur/src/tree_titer.py#L842) to try to catch common mistakes.

* Matching serum and virus strain names
	* Serum strain names were parsed from the first and second rows of the matrix, abbreviations for locations had to be added to `name_abbrev`.
	* Some strain names included extra information before the actual name (NYMCX-263BA/HK/4801/2014, IVR-159 (A/Victoria/502/2010), 1,3B/FLORIDA/4/2006). Used regex to match these but other patterns may arise when only the actual strain name is needed.
	* Some serum names included date information that needed to be reformatted (B/SHANDONG/JUL-97, A/CALIFORNIA/9-APR). Again used regex to reformat, but other patterns may arise.
	* `check_strain_names` looks for serum strain names that are potentially not parsed correctly by looking for a matching reference or test virus strain name. This helped me spot patterns of incorrectly parsed serum names.
* Titer measurement values
	* The HI assay uses two fold dilutions to measure antigenic similarity between viruses. So each titer measurement can only be certain values. I manually changed some values like '160' was most likely incorrectly entered and should be '180'.
	* The pdf converter sometimes made mistakes like combining two columns into one, half of one column into another (80 | 160 would become 8 | 0 160), regex handled some of these mistakes but others may arise, the function `check_titer_values` tries to spot these.
* Column labels
	* The HI tables have included more columns over time. They seems to always include `viruses`, `collection date` and `passage history`, but have added `genetic group` and `other information` columns.
	* `determine_columns` looks for field names from `self.table_column_names` in the first row of the HI table to label data parsed from that column. Some of the February 2016 column names were parsed to second row and this method didn't work in that case, need to manually move column names up in the csv files. Removed `other information` and blank columns from the HI matrix.

## Downloading

Measurements can be downloaded from tdb.

* Downloads all documents in database
* Prints result to designated `.tsv` or `.json` file.
	* Writes null attributes as '?'
	* Writes text file description in this order (0:`virus`, 1:`serum`, 2:`ferret_id`, 3:`source`, 4:`titer`)

### Commands

Command line arguments to run `download.py`:

* `-db --database`: database to download from, eg. `tdb`, `test_db`
* `-v --virus`: virus table to interact with, eg. `h1n1pdm`
* `--host`: host to be include in download, multiple arguments allowed
* `--path`: path to dump output files to, default is `data/`
* `--ftype`: output file format, default is `tsv`, other option is `json`
* `--fstem`: output file stem name, default is `VirusName\_Year\_Month\_Date`
* `--auth\_key`: authorization key for rethink database
* `--host`: rethink host url

### Examples

Download all H1N1pdm titers:

    python tdb/download.py -db tdb -v h1n1pdm

## Backup and Restore

TDB tables can be backed up to S3 or locally.

* Backups can be run manually or continuously everyday
* Backs up all tables in database
* Restoration keeps current documents in database but overwrites conflicting documents with the same primary key

### Examples

Backup `tdb.flu` to s3 backup file

	python tdb/backup.py -db tdb --backup_s3

Backup `tdb` to local backup file

	python tdb/backup.py -db tdb --backup_local

Backup `tdb` to s3 backup file everyday

	python tdb/backup.py -db tdb --continuous_backup --backup_s3

Restore `tdb.flu` from s3 backup file `2016-08-17_tdb_flu.tar.gz`

	python tdb/restore.py -db tdb -v flu --backup_s3 --restore_date 2016-08-17

## Append

Append documents to other tables in different databases. Useful for testing outcomes in a test database.table

### Examples

Append `tdb` flu documents to `test_tdb.flu`:

	python vdb/append.py -v flu --from_database tdb --to_database test_tdb


## Accessing the Database

All titer measurements are stored using [Rethinkdb deployed on AWS](https://www.rethinkdb.com/docs/paas/#deploying-on-aws). To access tdb you need an authorization key. This can be passed as a command line argument (see above) or set as an environment variable with a bash script, by running `source environment_rethink.sh`:

```shell
#!/bin/bash
export RETHINK_AUTH_KEY=EXAMPLE_KEY
export RETHINK_HOST=EXAMPLE_HOST
export NCBI_EMAIL=example\@email.org
```
