# TDB
The titer database (TDB) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database. 

## Uploading

## Schema

### Attribute Requirements

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