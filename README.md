# VDB

## Schema

* `Strain`: primary key. The canonical strain name. For flu this would be something like `A/Perth/16/2009`.
* `Virus`: Virus type. Loose term for like viruses (viruses that you'd want to include in a single tree). Examples include `H3N2`, `H1N1pdm`, `Vic`, `Yam`, `Ebola`, `Zika`.
* `Date`: Collection date in `YYYY-MM-DD` format, for example, `2016-02-28`.
* `Country`: Collection country in CamelCase format. See [here](https://github.com/blab/nextflu/blob/master/augur/source-data/geo_synonyms.tsv) for examples.
* `Division`: Administrative division in CamelCase format. Where available, Null otherwise.
* `Location`: Specific location in CamelCase format. Where available, Null otherwise.
* `Sequences`: list of...
  * `Accession`: Accession number. Where available, Null otherwise.
  * `Source`: Genbank, GISAID, etc...
  * `Locus`: gene or genomic region, `HA`, `NA`, `Genome`, etc...
  * `Sequence`: Actual sequence. Upper case.

