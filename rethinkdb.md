## Introduction to Rethinkdb

[Rethinkdb](https://www.rethinkdb.com/) stores items as JSON documents in a non-relational 
SQL ([noSQL](https://en.wikipedia.org/wiki/NoSQL)) type database. 
The noSQL database model is sometimes seen as more flexible than SQL databases where it's 
important to have a set schema from the beginning. Rethinkdb does not care how the 
data is formatted so data inserted later can have new information. Data cleaning and 
formatting must occur before inserting into the database if you want a certain schema enforced. 
Once it's in the correct format 
the [rethinkdb python module](https://www.rethinkdb.com/docs/guide/python/)
makes it very easy to insert and download JSON documents from 
the database as python dictionaries. 