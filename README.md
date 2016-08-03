## Introduction

The nextstrain project is an attempt to make flexible informatic pipelines and visualization tools to track ongoing pathogen evolution as sequence data emerges. The nextstrain project derives from [nextflu](https://github.com/blab/nextflu), which was specific to influenza evolution.

nextstrain is comprised of three components:

* [db](https://github.com/blab/nextstrain-db): database and IO scripts for sequence and serological data
* [augur](https://github.com/blab/nextstrain-augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/blab/nextstrain-auspice): web app to visualize resulting inferences

## db

The nextstrain db stores viral sequences and serological data in a [rethink database](RETHINKDB.md). The current database and scripts is designed around influenza and Zika viruses.

### vdb

The [virus database (vdb)](vdb/) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

### tdb

The [titer database (tdb)](tdb/) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database.

## Install

Clone the repo and load submodules:

    git clone https://github.com/blab/nextstrain-db.git
    git submodule update --init --recursive

Install Python bindings needed to run upload/download scripts:

    pip install rethinkdb==2.2.0.post2

Install Chateau Web UI:

    npm install

Backup and restore functionality requires the rethinkdb command line utility. This can be installed by following instructions [here](http://www.rethinkdb.com/docs/install/).

## Chateau

[Chateau](https://github.com/blab/chateau/) allows easy web access to the database. To run, do the following:

#### For remote rethink instance

1. Set environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY`.
2. Run with `npm run chateau` from directory `nextstrain-db/`.
3. Go to `http://localhost:3000/`.

#### For local rethink instance

2. Run with `npm run chateau-local` from directory `nextstrain-db/`.
3. Go to `http://localhost:3001/`.

Chateau configurations are stored in [`config.js`](config.js) for remote server and [`config_local.js`](config_local.js) for local server.

## License and copyright

Copyright 2016 Trevor Bedford and Charlton Callender.

Source code to nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
