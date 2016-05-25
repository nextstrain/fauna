## Introduction to Rethinkdb

[Rethinkdb](https://www.rethinkdb.com/) stores items as JSON documents in a non-relational 
SQL ([noSQL](https://en.wikipedia.org/wiki/NoSQL)) type database. This was useful for 
nextstrain because

* The noSQL database model is very flexible so data inserted later can have new 
information easy included. As a result though, more stringent data cleaning and formatting 
must occur before inserting into the database if you want a certain schema softly enforced. 
* Once it's in the correct format the [rethinkdb python module](https://www.rethinkdb.com/docs/guide/python/)
makes it very easy to insert and download JSON documents from the database as python 
dictionaries. 
* The opensource program [chateau](https://github.com/neumino/chateau) allows for easy
visualization and editing of the database. 

## Installation for nextstrain-db
1. Clone nextstrain-db repository `git clone https://github.com/blab/nextstrain-db.git` 
2. Install correct versions of rethinkdb and chateau from [package.json](package.json) `sudo npm install -g` 
3. Install correct version of rethinkdb python module `pip install rethinkdb==2.2.0.post2`

## Using rethinkdb

### Start a rethinkdb instance
The first thing you need to do to start using rethinkdb is [start a rethinkdb server](https://rethinkdb.com/docs/start-a-server/).
This can be done by just running `rethinkdb` in the terminal. The server can be hosted locally
or in the cloud (like on [AWS](https://rethinkdb.com/docs/paas/#deploying-on-aws)). If hosted
locally, data will be saved after stopping the server and loaded when rebooting the server.

### Import rethinkdb driver
```
import rethinkdb as r 
```

### Open connection to rethink database

### Create new databases and tables
Rethinkdb allows multiple databases within each instance and multiple tables within each
database. Each table stores similarly related data and each database stores similarly related
tables. It's recommended to have test databases and tables to see what the data looks like.

Chateau allows for easy creation of databases and tables. Sometimes it says it fails when
creating new tables but it actually does create them. When creating the table, you can specify 
what fields to use as the primary key. 

### Indexes

### Inserting documents into the database

### Downloading documents from the database

### Updating documents in the database