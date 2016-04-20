# TDB
The titer database (TDB) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database. 

## Uploading
Titer measurements from the Crick Worldwide Influenza Centre can be uploaded from csv files that are in a standard format.

## Schema

Each document in the database represents an HI test between a virus and serum from a specific ferret. Multiple titer measurements are stored if the test is repeated.

* `virus`: strain name of virus tested.
* `serum`: strain name serum was raised against. 
* `ferret_id`: id of the ferret the serum was raised in.
* `index`: Used as compound index in rethinkdb. Array with `virus`, `serum`, `ferret_id`.
* `titer`: List of all titer measurements for this test. 
* `source`: List of document names from which titer measurements were added. 
* `date_modified`:  Last modification date for document in `YYYY-MM-DD` format.
* `date`: Collection date of virus in `YYYY-MM-DD` format, for example, `2016-02-28`.
* `passage`: Passage history of the virus.
* `group`: Genetic group of the virus. 
* `ref`: Boolean for whether the virus is a reference virus. `True` if reference virus. 

### Attribute Requirements
Tests with null values for required attributes will be filtered out of those uploaded. Viruses with missing optional attributes will still be uploaded
* Required attributes: `virus`, `serum`, `ferret_id`, `index`, `titer`, `source`, `date_modified`
* Optional attributes: `date`, `passage`, `group`, `ref`

### Commands

### Examples:

## Downloading

### Commands

### Examples:

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