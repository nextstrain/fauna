import os, re, datetime, json, csv
from Bio import SeqIO
from Bio import Entrez
import requests

class parse(object):
    def __init__(self, **kwargs):
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
            viruses, sequences = self.get_entrez_viruses(gi, **kwargs)
        elif ftype is not None and fname is not None:
            if ftype == 'genbank':
                viruses, sequences = self.parse_gb_file(path + fname, **kwargs)
            elif ftype == 'accession':
                accessions = self.parse_accession_file(path + fname, **kwargs)
                self.entrez_email(email)
                gi = self.get_GIs(accessions)
                viruses, sequences = self.get_entrez_viruses(gi, **kwargs)
            elif ftype == 'fasta':
                viruses, sequences = self.parse_fasta_file(path + fname, **kwargs)
            elif ftype == 'tsv':
                viruses, sequences = self.parse_tsv_file(path + fname, **kwargs)
                print("Parsed " + str(len(viruses)) + " viruses and " + str(len(sequences)) + " sequences from file " + path+fname)
        else:
            raise Exception("No input file name and type defined or accessions given")
        return (viruses, sequences)

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

    def parse_fasta_file(self, fasta, virus_fasta_fields, sequence_fasta_fields, **kwargs):
        '''
        Parse FASTA file with default header formatting
        :return: list of documents(dictionaries of attributes) to upload
        '''
        header_fixes = False
        if (kwargs["fasta_header_fix"]):
            header_fixes = {}
            try:
                with open(kwargs["fasta_header_fix"], 'rU') as fh:
                    for line in fh:
                        if not line.startswith('#'):
                            k, v = line.strip().split("\t")
                            header_fixes[k] = v
            except IOError:
                raise Exception(kwargs["fasta_header_fix"], "not found")

        viruses = []
        sequences = []
        try:
            handle = open(fasta, 'r')
        except IOError:
            raise Exception(fasta, "not found")
        else:
            for record in SeqIO.parse(handle, "fasta"):
                if header_fixes:
                    try:
                        record.description = header_fixes[record.description]
                    except KeyError:
                        raise Exception(record.description, "not in header fix file. Fatal.")
                content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
                v = {key: content[ii] if ii < len(content) else "" for ii, key in virus_fasta_fields.items()}
                s = {key: content[ii] if ii < len(content) else "" for ii, key in sequence_fasta_fields.items()}
                s['sequence'] = str(record.seq).lower()
                v = self.add_virus_fields(v, **kwargs)
                s = self.add_sequence_fields(s, **kwargs)
                sequences.append(s)
                viruses.append(v)
            handle.close()
        return (viruses, sequences)

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
            raise Exception(tsv, "not found")
        else:
            with open(tsv) as infile:
                table_reader = csv.reader(infile, delimiter="\t")
                header = {i: element for i, element in enumerate(next(table_reader), 0)}
                for row in table_reader:
                    v = {key: row[ii] if ii < len(row) else "" for ii, key in header.items()}
                    self.add_virus_fields(v, **kwargs)
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
            raise Exception(gb, "not found")
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
            raise Exception(acc, "not found")
        else:
            accessions = []
            for acc in handle:
                accessions.append(acc.strip())
            return accessions

    def add_virus_fields(self, v, host, country, **kwargs):
        '''
        add fields to the viruses defined at the command line
        '''
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
        v['virus'] = self.virus
        v['timestamp'] = self.rethink_io.get_upload_timestamp()
        v['virus_inclusion_date'] = self.rethink_io.get_upload_date()
        v['sequences'] = []
        v['number_sequences'] = 0
        return v

    def add_sequence_fields(self, v, locus, authors, title, source, url, public=True, **kwargs):
        '''
        add fields to the sequences defined aat the commandline
        '''
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
        if 'title' not in v and title is not None:
            if title == 'null' or title == 'none' or title == 'None':
                v['title'] = None
            else:
                v['title'] = title
        if 'source' not in v and source is not None:
            if source == 'null' or source == 'none' or source == 'None':
                v['source'] = None
            else:
                v['source'] = source
        if 'url' not in v and url is not None:
            if url == 'null' or url == 'none' or url == 'None':
                v['url'] = None
            else:
                v['url'] = url
        if 'public' not in v and public is not None:
            v['public'] = public
        v['virus'] = self.virus
        v['timestamp'] = self.rethink_io.get_upload_timestamp()
        v['sequence_inclusion_date'] = self.rethink_io.get_upload_date()
        return v

    def get_GIs(self, accessions, n_entrez=2500, **kwargs):
        '''
        Use entrez esearch to get genbank identifiers from accession numbers
        '''
        retmax = 10**5  # max records to retrieve at once; 10^5 is documented limit, but >2500 reproducibly throws errors
        queries = []
        giList = []

        for i in sorted(range(0, len(accessions), n_entrez)): # split accessions list into 2500-long portions
            queries.append(" ".join(accessions[i:i+n_entrez])) # convert list to ' ' separated string

        assert sum([len(q.split()) for q in queries]) == len(accessions) # sanity check

        for q in queries:
            handle = Entrez.esearch(db=self.gbdb, term=q, retmax=retmax)    # retrieve xml of search results
            giList += Entrez.read(handle)['IdList'] # pull GI numbers from handle
        return giList

    def get_entrez_viruses(self, giList, **kwargs):
        '''
        Use entrez efetch to get genbank entries from genbank identifiers
        '''
        ## define batch size for download
        batchSize = 100

        # post NCBI query
        try:
            search_handle = Entrez.epost(db=self.gbdb, id=",".join(giList))
            search_results = Entrez.read(search_handle)
            webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
        except:
            print("ERROR: Couldn't connect with entrez, please run again")

        viruses = []
        sequences = []
        #fetch all results in batch of batchSize entries at once
        for start in range(0,len(giList),batchSize):
            #fetch entries in batch
            try:
                handle = Entrez.efetch(db=self.gbdb, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
            except IOError:
                print("ERROR: Couldn't connect with entrez, please run again")
            else:
                result = self.parse_gb_entries(handle, **kwargs)
                viruses.extend(result[0])
                sequences.extend(result[1])
        return (viruses, sequences)

    def parse_gb_entries(self, handle, **kwargs):
        '''
        Go through genbank records to get relevant virus information
        '''
        viruses, sequences = [], []
        SeqIO_records = SeqIO.parse(handle, "genbank")
        for record in SeqIO_records:
            v = {}
            s = {}
            s['source'] = 'genbank'
            s['accession'] = re.match(r'^([^.]*)', record.id).group(0).upper()  # get everything before the '.'?
            s['sequence'] = str(record.seq).lower()
            # set up as none and overwrite
            s["title"] = None
            s['authors'] = None
            s["puburl"] = None
            s["journal"] = None
            print("Processing genbank file for " + s['accession'])
            # all the references (i.e. the papers / direct-submission notes)
            references = record.annotations["references"]

            if len(references):
                # is there a reference which is not a "Direct Submission"?
                titles = [reference.title for reference in references]
                try:
                    idx = [i for i, j in enumerate(titles) if j is not None and j != "Direct Submission"][0]
                except IndexError: # fall back to direct submission
                    idx = [i for i, j in enumerate(titles) if j is not None][0]
                reference = references[idx] # <class 'Bio.SeqFeature.Reference'>
                keys = reference.__dict__.keys()
                s['title'] = reference.title
                if "authors" in keys and reference.authors:
                    first_author = re.match(r'^([^,]*)', reference.authors).group(0)
                    s['authors'] = first_author + " et al"
                if "journal" in keys and reference.journal:
                    s['journal'] = reference.journal
                if "pubmed_id" in keys and reference.pubmed_id:
                    s["puburl"] = "https://www.ncbi.nlm.nih.gov/pubmed/" + reference.pubmed_id
            else:
                print("Couldn't find the reference for " + s['accession'])

            # print(" *** Accession: {} title: {} authors: {} journal: {} paperURL: {}".format(s['accession'], s['title'], s['authors'], s['journal'], s['puburl']))

            s['url'] = "https://www.ncbi.nlm.nih.gov/nuccore/" + s['accession']
            #s['url'] = self.get_doi_url(url, s['title'], first_author)

            record_features = record.features
            for feat in record_features:
                if feat.type == 'source':
                    qualifiers = feat.qualifiers
                    if 'collection_date' in qualifiers:
                        v['collection_date'] = self.convert_gb_date(qualifiers['collection_date'][0])
                    if 'country' in qualifiers:
                        v['country'] = re.match(r'^([^:]*)', qualifiers['country'][0]).group(0)
                    if 'strain' in qualifiers:
                        v['strain'] = qualifiers['strain'][0]
                        s['strain'] = qualifiers['strain'][0]
                    elif 'isolate' in qualifiers:
                        v['strain'] = qualifiers['isolate'][0]
                        s['strain'] = qualifiers['isolate'][0]
                    else:
                        print("Couldn't parse strain name for " + s['accession'])
            v = self.add_virus_fields(v, **kwargs)
            s = self.add_sequence_fields(s, **kwargs)
            viruses.append(v)
            sequences.append(s)
            #self.fix_casing(v)
            #self.fix_boolean(v)
        handle.close()
        print(str(len(viruses)) + " genbank entries parsed")
        return (viruses, sequences)

    def get_doi_url(self, url, title, author):
        '''
        Use crossref api to look for matching title and author name to link to DOI
        Depcrecated for the moment. Crossref was hanging randomly.
        '''
        if title is not None:
            num = str(2)
            print("Calling crossref for " + url)
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
            if collection_date.split('-')[0].isdigit():
                if int(collection_date.split('-')[0]) < 13:
                    return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%m-%Y'), '%Y-%m-XX')
                else:
                    return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%Y-%m'), '%Y-%m-XX')
            else:
                return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%b-%Y'), '%Y-%m-XX')
        elif N_fields == 3:
            if int(collection_date.split('-')[0]) < 32:
                return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%d-%b-%Y'), '%Y-%m-%d')
            else:
                return datetime.datetime.strftime(datetime.datetime.strptime(collection_date,'%Y-%m-%d'), '%Y-%m-%d')

    def get_upload_date(self):
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))
