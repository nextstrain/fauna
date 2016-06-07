import os, re, datetime, json
from Bio import SeqIO
from Bio import Entrez
import requests

class parse(object):
    def __init__(self, **kwargs):
        if 'fasta_fields' in kwargs:
            self.fasta_fields = kwargs['fasta_fields']
        if 'accessions' in kwargs:
            self.accessions = kwargs['accessions']
        self.gbdb = "nuccore"

    def entrez_email(self, email, **kwargs):
        #define email for entrez login
        if email is not None:
            self.email = email
        elif 'NCBI_EMAIL' in os.environ:
            self.email = os.environ['NCBI_EMAIL']
        else:
            raise Exception("Missing NCBI email")
        Entrez.email = self.email

    def parse(self, path, fname, ftype, email, **kwargs):
        if self.accessions is not None:
            accessions = [acc.strip() for acc in self.accessions.split(",")]
            self.entrez_email(email)
            gi = self.get_GIs(accessions)
            self.viruses = self.get_entrez_viruses(gi, **kwargs)
        elif ftype is not None and fname is not None:
            if ftype == 'genbank':
                self.viruses = self.parse_gb_file(path + fname, **kwargs)
            elif ftype == 'accession':
                accessions = self.parse_accession_file(path + fname, **kwargs)
                self.entrez_email(email)
                gi = self.get_GIs(accessions)
                self.viruses = self.get_entrez_viruses(gi **kwargs)
            elif ftype == 'fasta':
                self.viruses = self.parse_fasta_file(path + fname, **kwargs)
            elif ftype == 'tsv':
                self.viruses = self.parse_tsv_file(path + fname, **kwargs)                
        else:
            raise Exception("No input file name and type defined or accessions given")

    def add_other_attributes(self, v, locus, authors, host, country, source, public=True, **kwargs):
        '''
        Add attributes to all viruses to be uploaded that are included at the command line
        '''
        if 'strain' in v:
            v['strain'] = self.fix_name(v['strain'])
        for field in self.grouping_upload_fields + self.grouping_optional_fields:
            if field in v:
                pass
            elif field in kwargs and kwargs[field] is not None:
                v[field] = kwargs[field]
        v['virus'] = self.virus
        v['timestamp'] = self.rethink_io.get_upload_timestamp()
        if 'locus' not in v and locus is not None:
            if locus == 'null' or locus == 'none' or locus == 'None':
                v['locus'] = None
            else:
                v['locus'] = locus
        if 'authors' not in v and authors is not None:
            if authors == 'null' or authors == 'none' or authors == 'None':
                v['authors'] = None
            else:
                v['authors'] = authors
        if 'host' not in v and host is not None:
            if host == 'null' or host == 'none' or host == 'None':
                v['host'] = None
            else:
                v['host'] = host
        if 'country' not in v and country is not None:
            if country == 'null' or country == 'none' or country == 'None':
                v['country'] = None
            else :
                v['country'] = country
        if 'source' not in v and source is not None:
            if source == 'null' or source == 'none' or source == 'None':
                v['source'] = None
            else:
                v['source'] = source
        if 'public' not in v and public is not None:
            v['public'] = public

    def fix_casing(self, v):
        '''
        force lower case on fields besides strain, title, authors
        '''
        for field in v:
            if field in ['strain', 'title', 'authors', 'accession']:
                pass
            elif v[field] is not None and isinstance(v[field], str):
                v[field] = v[field] = v[field].lower().replace(' ', '_')

    def fix_boolean(self, v):
        '''
        replace strings of "true" and "false" with proper booleans
        '''
        for field in v:
            if v[field] is not None:
                if v[field] == "true":
                    v[field] = True
                if v[field] == "false":
                    v[field] = False

    def parse_fasta_file(self, fasta, **kwargs):
        '''
        Parse FASTA file with default header formatting
        :return: list of documents(dictionaries of attributes) to upload
        '''
        viruses = []
        try:
            handle = open(fasta, 'r')
        except IOError:
            print(fasta, "not found")
        else:
            for record in SeqIO.parse(handle, "fasta"):
                content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
                v = {key: content[ii] if ii < len(content) else "" for ii, key in self.fasta_fields.items()}
                v['sequence'] = str(record.seq)
                self.add_other_attributes(v, **kwargs)
                self.fix_casing(v)
                self.fix_boolean(v)                
                viruses.append(v)
            handle.close()
        return viruses

    def parse_tsv_file(self, tsv, **kwargs):
        '''
        Parse TSV file with default field ordering
        This is the same ordering as specified in 'fasta_fields'
        :return: list of documents(dictionaries of attributes) to upload
        '''
        import csv
        viruses = []
        try:
            os.path.isfile(tsv)
        except IOError:
            print(tsv, "not found")
        else:
            with open(tsv) as infile:
                table_reader = csv.reader(infile, delimiter="\t")
                header = {i: element for i, element in enumerate(next(table_reader), 0)}
                for row in table_reader:
                    v = {key: row[ii] if ii < len(row) else "" for ii, key in header.items()}
                    self.add_other_attributes(v, **kwargs)
                    self.fix_casing(v)
                    self.fix_boolean(v)
                    viruses.append(v)
        return viruses

    def parse_gb_file(self, gb, **kwargs):
        '''
        Parse genbank file
        :return: list of documents(dictionaries of attributes) to upload
        '''
        try:
            handle = open(gb, 'r')
        except IOError:
            print(gb, "not found")
        else:
            return self.parse_gb_entries(handle, **kwargs)

    def parse_accession_file(self, acc, **kwargs):
        '''
        Parse file for list of accession numbers to be uploaded to vdb
        :return: list of documents(dictionaries of attributes) to upload
        '''
        try:
            handle = open(acc, 'r')
        except IOError:
            print(acc, "not found")
        else:
            accessions = []
            for acc in handle:
                accessions.append(acc.strip())
            return accessions

    def get_GIs(self, accessions, **kwargs):
        '''
        Use entrez esearch to get genbank identifiers from accession numbers
        '''
        retmax = 10**9
        query = " ".join(accessions)
        handle = Entrez.esearch(db=self.gbdb, term=query, retmax=retmax)
        giList = Entrez.read(handle)['IdList']
        return giList

    def get_entrez_viruses(self, giList, **kwargs):
        '''
        Use entrez efetch to get genbank entries from genbank identifiers
        '''
        ## define batch size for download
        batchSize = 100

        # post NCBI query
        search_handle = Entrez.epost(db=self.gbdb, id=",".join(giList))
        search_results = Entrez.read(search_handle)
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]

        viruses = []
        #fetch all results in batch of batchSize entries at once
        for start in range(0,len(giList),batchSize):
            #fetch entries in batch
            try:
                handle = Entrez.efetch(db=self.gbdb, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
            except IOError:
                print("Couldn't connect with entrez")
            else:
                viruses.extend(self.parse_gb_entries(handle, **kwargs))
        return viruses

    def parse_gb_entries(self, handle, **kwargs):
        '''
        Go through genbank records to get relevant virus information
        '''
        viruses = []
        for record in SeqIO.parse(handle, "genbank"):
            v = {}
            v['source'] = 'genbank'
            v['accession'] = re.match(r'^([^.]*)', record.id).group(0).upper()  # get everything before the '.'?
            v['sequence'] = str(record.seq)
            reference = record.annotations["references"][0]
            if reference.title is not None and reference.title != "Direct Submission":
                v['title'] = reference.title
            else:
                print("Couldn't find reference title for " + v['accession'])
                v['title'] = None
            if reference.authors is not None:
                first_author = re.match(r'^([^,]*)', reference.authors).group(0)
            else:
                print("Couldn't parse authors for " + v['accession'])
                first_author = None
            url = "http://www.ncbi.nlm.nih.gov/nuccore/" + v['accession']
            v['url'] = self.get_gb_url(url, v['title'], first_author)
            v['authors'] = first_author + " et al"

            record_features = record.features
            for feat in record_features:
                if feat.type == 'source':
                    qualifiers = feat.qualifiers
                    v['date'] = self.convert_gb_date(qualifiers['collection_date'][0])
                    v['country'] = re.match(r'^([^:]*)', qualifiers['country'][0]).group(0)
                    if 'strain' in qualifiers:
                        v['strain'] = qualifiers['strain'][0]
                    elif 'isolate' in qualifiers:
                        v['strain'] = qualifiers['isolate'][0]
                    else:
                        print("Couldn't parse strain name for " + v['accession'])
            self.add_other_attributes(v, **kwargs)
            self.fix_casing(v)
            self.fix_boolean(v)
            viruses.append(v)
        handle.close()
        print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return viruses

    def get_gb_url(self, url, title, author):
        '''
        Use crossref api to look for matching title and author name to link to DOI"
        '''
        if title is not None:
            num = str(2)
            response = json.loads(requests.get('http://api.crossref.org/works?query=%' + title + '%22&rows=' + num).text)
            items = response['message']['items']
            for item in items:
                if 'title' in item and item['title'][0] == title:
                    if 'author' in item and item['author'][0]['family'] == author:
                        if 'DOI' in item:
                            url = 'http://dx.doi.org/' + item['DOI']
        return url

    def convert_gb_date(self, collection_date):
        '''
        Converts calendar dates between given formats
        '''
        N_fields = len(collection_date.split('-'))
        if N_fields == 1:
            return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%Y'), '%Y-XX-XX')
        elif N_fields == 2:
            return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%b-%Y'), '%Y-%m-XX')
        elif N_fields == 3:
            return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%d-%b-%Y'), '%Y-%m-%d')

    def get_upload_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))
