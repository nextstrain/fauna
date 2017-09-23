## Introduction

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data. We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread. Our goal is to aid epidemiological understanding and improve outbreak response.

Nextstrain is comprised of three primary components:

* [fauna](https://github.com/nextstrain/fauna): database and IO scripts for sequence and serological data
* [augur](https://github.com/nextstrain/augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/nextstrain/auspice): web app to visualize resulting inferences

Resulting data and inferences are available live at the website [nextstrain.org](http://nextstrain.org).

## fauna

*Definition: The animals of a given region or period considered as a whole. Also, prophetic Roman deity.*

The fauna database stores viral sequences and serological data in [RethinkDB](RETHINKDB.md). The current database and scripts are designed around influenza, Ebola and Zika viruses, but with the intention of provided a general purpose tool.

_Note: In most cases, it will be easier to pass augur a self-prepared FASTA file than to use fauna with the overhead of launching a RethinkDB instance. If you are new to Nextstrain, we suggest you skip fauna and proceed to directly to augur._

### vdb

The [virus database (vdb)](vdb/) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

### tdb

The [titer database (tdb)](tdb/) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database.

### Supported virus builds

We maintain notes on [supported virus builds](builds/).

## Install

Clone the repo and load submodules:

    git clone https://github.com/nextstrain/fauna.git
    git submodule update --init --recursive

Install Python modules needed to run upload/download scripts:

    pip install -r requirements.txt

Install Chateau Web UI:

    npm install

Backup and restore functionality requires the rethinkdb command line utility. This can be installed by following instructions [here](http://www.rethinkdb.com/docs/install/). With Homebrew, you can just do:

    brew install rethinkdb

## Chateau

[Chateau](https://github.com/nextstrain/chateau/) allows easy web access to the database. To run, do the following:

#### For remote rethink instance

1. Set environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY`.
2. Run with `npm run chateau` from directory `fauna/`.
3. Go to `http://localhost:3000/`.

#### For local rethink instance

2. Run with `npm run chateau-local` from directory `fauna/`.
3. Go to `http://localhost:3001/`.

Chateau configurations are stored in [`config.js`](config.js) for remote server and [`config_local.js`](config_local.js) for local server.

## License and copyright

Copyright 2016-2017 Trevor Bedford.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
