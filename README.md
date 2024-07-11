## Introduction

Nextstrain is an open-source project to harness the scientific and public health potential of pathogen genome data. We provide a continually-updated view of publicly available data with powerful analytics and visualizations showing pathogen evolution and epidemic spread. Our goal is to aid epidemiological understanding and improve outbreak response.

Resulting data and inferences are available live at the website [nextstrain.org](https://nextstrain.org). Documentation is available at [nextstrain.org/docs](https://nextstrain.org/docs).

## Fauna

*Definition: The animals of a given region or period considered as a whole. Also, prophetic Roman deity.*

The Fauna database stores viral sequences and serological data in [RethinkDB](RETHINKDB.md). The current database and scripts are designed around influenza, Ebola and Zika viruses, but with the intention of provided a general purpose tool.

_Note: In most cases, it will be easier to pass augur a self-prepared FASTA file than to use Fauna with the overhead of launching a RethinkDB instance. If you are new to Nextstrain, we suggest you skip Fauna and proceed to directly to Augur._

_Note: The Nextstrain team's own internal Fauna instance is private because it contains data that cannot be public.  If you're part of the Nextstrain team, you can see [details of its current setup on AWS](https://github.com/nextstrain/private/issues/86#issuecomment-1793247244)._

### vdb

The [virus database (vdb)](vdb/) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

### tdb

The [titer database (tdb)](tdb/) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database.

### Supported virus builds

We maintain notes on [supported virus builds](builds/).

## Install

Fauna requires Python 3.
> **_note_**: A subset of tdb upload scripts use `xlrd` (v1.2.0) which will not work with Python 3.8 or newer

Clone the repo and load submodules:

    git clone https://github.com/nextstrain/fauna.git
    cd fauna
    git submodule update --init --recursive

Install Python modules needed to run upload/download scripts:

    python3 -m pip install -r requirements.txt

Install Chateau Web UI:
> **_note_**: this step is optional and only required if you plan to use Chateau to explore the data

    cd chateau
    npm install --production

Backup and restore functionality requires the rethinkdb command line utility. This can be installed by following instructions [here](http://www.rethinkdb.com/docs/install/). With Homebrew, you can just do:

    cd ..
    brew install rethinkdb

Most functions have been converted to work only in Python 3. However particular calls may not have been converted. If you run into a Python 2/3 error please note in an issue.

## Environment variables

Throughout we assume the existence of environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY`. We do not share these variables here, because for security reasons our RethinkDB instance is private. To use Fauna you will need to set up your own RethinkDB instance as described [here](RETHINKDB.md). This instance can be local, in which case variables will be:

* `RETHINK_HOST`: `localhost`
* `RETHINK_AUTH_KEY`: ``

Or this instance can be remote, in which case follow the [RethinkDB docs to configure](https://rethinkdb.com/docs/security/#securing-the-web-interface). Note that admin password is synonymous with RethinkDB `auth_key`.

## Chateau

[Chateau](https://github.com/neumino/chateau) (forked to [nextstrain/chateau/](https://github.com/nextstrain/chateau/tree/trs/tls-support)) allows easy web access to the database. To run, do the following:

1. Start a local rethinkdb server by running `rethinkdb`, then switch to a new terminal.
2. Set environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY`. If running locally set `RETHINK_HOST` to `localhost` and leave `RETHINK_AUTH_KEY` empty.
3. Run with `npm run chateau` from directory `fauna/`.
4. Go to `http://localhost:3000/`.

Chateau configurations are stored in [`config.js`](config.js).

## License and copyright

Copyright 2016-2020 Trevor Bedford.

Source code to Nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). Nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
