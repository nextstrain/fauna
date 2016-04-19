import os, re, datetime, json
from Bio import SeqIO
from Bio import Entrez
import requests

class vdb_parse(object):
    def __init__(self, **kwargs):
        if 'fasta_fields' in kwargs:
            self.fasta_fields = kwargs['fasta_fields']
        if 'accessions' in kwargs:
            self.accessions = kwargs['accessions']
        self.gbdb = "nuccore"

        #define email for entrez login
        if 'email' in kwargs:
            self.email = kwargs['email']
        if 'NCBI_EMAIL' in os.environ and self.email is None:
            self.email = os.environ['NCBI_EMAIL']
        if self.email is None:
            raise Exception("Missing NCBI email")
        Entrez.email = self.email

    def parse(self):
        if self.accessions is not None:
            accessions = [acc.strip() for acc in self.accessions.split(",")]
            gi = self.get_GIs(accessions)
            self.viruses = self.get_entrez_viruses(gi)
        elif self.ftype is not None and self.fname is not None:
            if self.ftype == 'genbank':
                self.viruses = self.parse_gb_file(self.path + self.fname)
            elif self.ftype == 'accession':
                accessions = self.parse_accession_file(self.path + self.fname)
                gi = self.get_GIs(accessions)
                self.viruses = self.get_entrez_viruses(gi)
            elif self.ftype == 'fasta':
                self.viruses = self.parse_fasta_file(self.path + self.fname)
        elif self.auto_upload:
            gi = self.auto_gb_upload()
            print(len(gi))
            raise Exception("stop")
            self.viruses = self.get_entrez_viruses(gi)

        else:
            raise Exception("No input file name and type defined or accessions given")

    def auto_gb_upload(self):
        '''

        :return:
        '''
        retmax = 10**9
        organism = "\"Zika virus\"[porgn]"
        start_date = "2015/01/01"
        end_date = self.get_upload_date()
        modification_date = "(\"" + start_date + "\"[MDAT] : \"" + end_date + "\"[MDAT])"
        minimum_length = 10000
        sequence_length = "(\"" + str(minimum_length) + "\"[SLEN] : \"" + str(minimum_length**2) + "\"[SLEN])"
        query = " AND ".join([organism, modification_date, sequence_length])
        handle = Entrez.esearch(db=self.gbdb, term=query, retmax=retmax)
        giList = Entrez.read(handle)['IdList']
        return giList

    def parse_fasta_file(self, fasta):
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
                v['sequence'] = str(record.seq).upper()
                v['virus'] = self.virus
                v['date_modified'] = self.get_upload_date()
                if 'locus' not in v and self.locus is not None:
                    v['locus'] = self.locus.title()
                if 'authors' not in v and self.authors is not None:
                    v['authors'] = self.authors.title()
                if 'subtype' not in v and self.vsubtype is not None:
                    v['subtype'] = self.vsubtype.title()
                if 'source' not in v and self.virus_source is not None:
                    v['source'] = self.virus_source.title()

                viruses.append(v)
            handle.close()
            print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return viruses

    def parse_gb_file(self, gb):
        '''
        Parse genbank file
        :return: list of documents(dictionaries of attributes) to upload
        '''
        try:
            handle = open(gb, 'r')
        except IOError:
            print(gb, "not found")
        else:
            return self.parse_gb_entries(handle)

    def parse_accession_file(self, acc):
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

    def get_GIs(self, accessions):
        '''
        Use entrez esearch to get genbank identifiers from accession numbers
        '''
        retmax = 10**9

        query = " ".join(accessions)
        handle = Entrez.esearch(db=self.gbdb, term=query, retmax=retmax)
        giList = Entrez.read(handle)['IdList']
        return giList

    def get_entrez_viruses(self, giList):
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
                viruses.extend(self.parse_gb_entries(handle))
        return viruses

    def parse_gb_entries(self, handle):
        '''
        Go through genbank records to get relevant virus information
        '''
        viruses = []
        for record in SeqIO.parse(handle, "genbank"):
            v = {}
            v['source'] = 'Genbank'
            v['accession'] = re.match(r'^([^.]*)', record.id).group(0).upper()  # get everything before the '.'?
            v['sequence'] = str(record.seq).upper()
            v['virus'] = self.virus
            v['date_modified'] = self.get_upload_date()
            reference = record.annotations["references"][0]
            if reference.title is not None and reference.title != "Direct Submission":
                v['title'] = reference.title
            else:
                print("Couldn't find reference title for " + v['accession'])
                v['title'] = None
            if reference.authors is not None:
                first_author = re.match(r'^([^,]*)', reference.authors).group(0).title()
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
                    if 'isolate' in qualifiers:
                        v['strain'] = qualifiers['isolate'][0]
                    elif 'strain' in qualifiers:
                        v['strain'] = qualifiers['strain'][0]
                    else:
                        print("Couldn't parse strain name for " + v['accession'])
            if 'locus' not in v and self.locus is not None:
                v['locus'] = self.locus.title()
            if 'authors' not in v and self.authors is not None:
                v['authors'] = self.authors.title()
            if 'subtype' not in v and self.vsubtype is not None:
                v['subtype'] = self.vsubtype.title()
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
