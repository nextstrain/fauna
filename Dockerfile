# nextstrain-db dockerfile

FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv
RUN apt-get install -y python-numpy python-scipy 
RUN apt-get install -y libpng-dev libfreetype6-dev pkg-config 

# python dependencies
RUN pip install rethinkdb==2.2.0.post2
RUN pip install biopython==1.64
RUN pip install geopy==1.11.0

# nextstrain-db
RUN git clone https://github.com/blab/nextstrain-db.git /db
WORKDIR /db/

# update
ADD http://www.timeapi.org/utc/now /tmp/bustcache
RUN git pull

# default process
CMD python -u vdb/flu_download.py -db test_vdb -v flu --select locus:HA subtype:h3n2 --interval collection_date:2016-01-01,2016-01-15 --pick_longest --fstem h3n2
