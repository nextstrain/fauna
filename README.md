# VDB

## Schema

* `Strain`: primary key. The canonical strain name. For flu this would be something like `A/Perth/16/2009`.
* `Date`: Collection date in `YYYY-MM-DD` format, for example, `2016-02-28`.
* `Country`: Collection country in CamelCase format. See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_synonyms.tsv) for examples.
* `Division`: Administrative division in CamelCase format. Where available, Null otherwise.
* `Location`: Specific location in CamelCase format. Where available, Null otherwise.
* `Sequences`: list of...
  * `Accession`: Accession number. Where available, Null otherwise.
  * `Database`: Genbank, GISAID, etc...
  * `Region`: gene or genomic region, `HA`, `NA`, `Genome`, etc...
  * `Sequence`: Actual sequence. Upper case.

