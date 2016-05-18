## Introduction

The nextstrain project is an attempt to make flexible informatic pipelines and visualization tools to track ongoing pathogen evolution as sequence data emerges. The nextstrain project derives from [nextflu](https://github.com/blab/nextflu), which was specific to influenza evolution.

nextstrain is comprised of three components:

* [db](https://github.com/blab/nextstrain-db): database and IO scripts for sequence and serological data
* [augur](https://github.com/blab/nextstrain-augur): informatic pipelines to conduct inferences from raw data
* [auspice](https://github.com/blab/nextstrain-auspice): web app to visualize resulting inferences

## Install

To install Python bindings, run:

    pip install rethinkdb==2.2.0.post2

To install Chateau Web UI, from the base directory of `nextstrain-db/`, run:

    sudo npm install -g

## db

The nextstrain db stores viral sequences and serological data for (currently) influenza and Zika viruses.

### vdb

The [virus database (vdb)](vdb/) is used to store viral information in an organized schema. This allows easy storage and querying of viruses which can be downloaded in formatted fasta or json files.

### tdb

The [titer database (tdb)](tdb/) is used to store titer measurements in an organized schema. This allows easy storage and downloading of all measurements in the database.

## Chateau

[Chateau](https://github.com/neumino/chateau/) allows easy web access to the database. To run, do the following:

1. Navigate to the base directory of `nextstrain-db/`.
2. Set environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY`.
3. Run with `chateau`.
4. Fire up a browser and go to `http://localhost:3000/`.

Chateau configuration is stored in the file [`config.js`](config.js).

## License and copyright

Copyright 2016 Trevor Bedford and Charlton Callender.

Source code to nextstrain is made available under the terms of the [GNU Affero General Public License](LICENSE.txt) (AGPL). nextstrain is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
