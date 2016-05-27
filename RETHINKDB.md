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
2. Install correct versions of rethinkdb and chateau from [package.json](package.json) with `sudo npm install -g` 
3. Install correct version of rethinkdb python module `pip install rethinkdb==2.2.0.post2`

## Using rethinkdb

### Start a rethinkdb instance
The first thing you need to do to start using rethinkdb is [start a rethinkdb server](https://rethinkdb.com/docs/start-a-server/).
This can be done by just running `rethinkdb` in the terminal. The server can be hosted locally
or in the cloud (like on [AWS](https://rethinkdb.com/docs/paas/#deploying-on-aws)). If hosted
locally, data will be saved after stopping the server and loaded when rebooting the server, but you must run `rehtinkdb` from the same directory.

### Import rethinkdb driver
```
import rethinkdb as r 
```

### Open connection to rethink database
To open a connection to a database on a local host run
```
r.connect(host='localhost', port=28015, db=database).repl()
```
To open a connection to a database on an external host that requires an authorization key run
```
r.connect(host=rethink_host, port=28015, db=database, auth_key=auth_key).repl()
```

### Create new databases and tables
Rethinkdb allows multiple databases within each instance and multiple tables within each
database. Each table stores similarly related data and each database stores similarly related
tables. It's recommended to have test databases and tables to see what the data looks like.

Chateau allows for easy creation of databases and tables. Sometimes it says it fails when
creating new tables but it actually does create them. When creating the table, you can specify 
what fields to use as the primary key. 

### Indexes
Each document in a table must have a unique primary key. This primary key can be a simple 
index (a string, in [vdb](vdb) this is `strain`) or a compound index (list of strings, 
used in [tdb](tdb)). If a primary key is not defined for a table, rethink will automatically 
assign a unique id to each document. 

### Inserting documents into the database
To [insert](https://rethinkdb.com/api/python/insert/) one JSON document or a whole list 
of documents into a specific table, run 
```
r.table(table).insert(documents).run()
```
By default, this will error if a document with the same primary key is already in the table.
There is also the option to 'replace' the old document, or 'update' fields of the old document
with the new fields. 
```
r.table(table).insert(documents, conflict='replace').run()
r.table(table).insert(documents, conflict='update').run()
```
### Downloading documents from the database
To download a cursor of all documents in a table you can run the following and loop 
through the iterable cursor
```
cursor = r.table(table).run()
for document in cursor:
    print(document)
```
It's also possible to filter the values you receive through rethinkdb or this can be done
after receiving all the documents in the table. 
```
cursor = r.table(table).filter(r.row[field] == field_value).run()
```
You can also retrieve a specific document using it's index. This returns `None` if the document
is not in the table.
```
r.db(database).table(table).get(index).run()
```

### More Information about Rethinkdb
[Ten-minute guide with RethinkDB and Python](https://www.rethinkdb.com/docs/guide/python/)

[Python ReQL command reference](https://www.rethinkdb.com/api/python/)