# VDB
The virus database (VDB) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

## Schema

* `Strain`: primary key. The canonical strain name. For flu this would be something like `A/Perth/16/2009`.
* `Virus`: Virus type in CamelCase format. Loose term for like viruses (viruses that you'd want to include in a single tree). Examples include `Flu`, `Ebola`, `Zika`.
* `Subtype`: Virus subtype in CamelCase format, where available, Null otherwise. `H3N2`, `H1N1pdm`, `Vic`, `Yam`
* `Date`: Collection date in `YYYY-MM-DD` format, for example, `2016-02-28`.
* `Region`: Collection region in CamelCase format.  See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_regions.tsv) for examples. 
* `Country`: Collection country in CamelCase format. See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_synonyms.tsv) for examples.
* `Division`: Administrative division in CamelCase format. Where available, Null otherwise.
* `Location`: Specific location in CamelCase format. Where available, Null otherwise.
* `Sequences`: list of...
  * `Accession`: Accession number. Where available, Null otherwise.
  * `Source`: Genbank, GISAID, etc... in CamelCase format.
  * `Authors`: Authors to attribute credit to. Where available, Null otherwise. in CamelCase format.
  * `Locus`: gene or genomic region, `HA`, `NA`, `Genome`, etc... in CamelCase format.
  * `Sequence`: Actual sequence. Upper case.

## Accessing the Database
All viruses are stored using [Rethinkdb deployed on AWS](https://www.rethinkdb.com/docs/paas/#deploying-on-aws)

To access vdb you need an authorization key. This can be passed as a command line argument (see below) or set as an environment variable with a bash script.

`source environment_rethink.sh`
```shell
#!/bin/bash
export RETHINK\_AUTH\_KEY=EXAMPLE\_KEY
```

## Uploading
Sequences can be uploaded from a fasta file to a virus specific table within vdb. It currently
* Reads a fasta file for description and sequence information.
	* Reads fasta description in this order for Zika\_vdb\_upload (0:`accession`, 1:`strain`, 2:`date`, 4:`country`, 5:`division`, 6:`location`)
	* Reading order can easily be changed in source code for different viruses or files
* Formats information to fit with vdb schema. 
* Uploads information to virus table
	* Appends to list of sequences if new accession number. If no accession number, appends if new sequence.
	* If strain already in database, default is to update attributes with new information only if current attribute is null
	
### Attribute Requirements
Viruses with null values for required attributes will be filtered out of those uploaded. Viruses with missing optional attributes will still be uploaded
* Required virus attributes: `strain`, `date`, `country`, `sequences`, `virus`, 
* Required sequence attributes: `source`, `locus`, `sequence`
* Optional virus attributes: `subtype`, `division`, `location`
* Optional sequence attributes: `accession`, `authors`
(`subtype` is required for Flu)

### Commands
Command line arguments to run vdb_upload:
* -db --database default='vdb', help=database to upload to. Ex 'vdb', 'test'
* -v --virus help=virus table to interact with. Ex 'Zika', 'Flu'
* --fname help=input fasta file name
* --source
* --locus
* --authors
* --subtype
* --overwrite default=False help=whether to overwrite existing non-null fields
* --path help=path to fasta file, default is data/virus/
* --auth\_key help=authorization key for rethink database

Example commands:

`python Flu_vdb_upload.py -db test -v flu --fname H3N2_gisaid_epiflu_sequence.fasta --source gisaid --subtype H3N2`

`python Zika_vdb_upload.py --database vdb --virus zika --fname zika_virological_02-22-2016.fasta --source Virological --locus Genome`

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
* -v --virus help=virus table to interact with. Ex 'Zika', 'Flu'
* --path help=path to dump output files to, default is data/
* --ftype help=output file format, default is 'fasta', other option is 'json'
* --fstem help=output file stem name, default is VirusName\_Year\_Month\_Date
* --auth\_key help=authorization key for rethink database

Example commands:

`python vdb_download.py -db vdb -v Zika --fstem Zika_test`
