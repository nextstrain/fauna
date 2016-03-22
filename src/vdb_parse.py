import os, re, datetime
from Bio import SeqIO
from Bio import Entrez

class vdb_parse(object):
    def __init__(self, fasta_fields, **kwargs):
        self.fasta_fields = fasta_fields
        self.accessions = kwargs['accessions']

    def parse(self):
        if self.ftype is not None and self.fname is not None:
            if self.ftype == 'genbank':
                self.viruses = self.parse_gb_file(self.path + self.fname)
            elif self.ftype == 'accession':
                accessions = self.parse_accession_file(self.path + self.fname)
                self.viruses = self.access_entrez(accessions)
            elif self.ftype == 'fasta':
                self.viruses = self.parse_fasta_file(self.path + self.fname)
        elif self.accessions is not None:
            accessions = [acc.strip() for acc in self.accessions.split(",")]
            self.viruses = self.access_entrez(accessions)
        else:
            raise Exception("No input file name and type defined or accessions given")

        return self.viruses

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

    def access_entrez(self, accessions):
        '''
        Use entrez database to get genbank entries for list of accessions
        '''
        #define email for entrez login
        db = "nuccore"
        if self.email is not None:
            Entrez.email = self.email
        else:
            raise Exception("Need to give email to upload via accession number")

        ## define batch size for download
        batchSize = 100
        retmax = 10**9

        # Get genbank identifiers from accession numbers using esearch
        query = " ".join(accessions)
        handle = Entrez.esearch(db=db, term=query, retmax=retmax)
        giList = Entrez.read(handle)['IdList']

        # post NCBI query
        search_handle = Entrez.epost(db=db, id=",".join(giList))
        search_results = Entrez.read(search_handle)
        webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]

        viruses = []
        #fetch all results in batch of batchSize entries at once
        for start in range(0,len(giList),batchSize):
            #fetch entries in batch
            try:
                handle = Entrez.efetch(db=db, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
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
            if reference.authors is not None:
                first_author = re.match(r'^([^,]*)', reference.authors).group(0).title()
                v['authors'] = first_author + " et al"
            else:
                print("Couldn't parse authors for " + v['accession'])
            if reference.title is not None:
                v['title'] = reference.title
            else:
                print("Couldn't parse reference title for " + v['accession'])
            v['url'] = "http://www.ncbi.nlm.nih.gov/nuccore/" + v['accession']

            record_features = record.features
            for feat in record_features:
                if feat.type == 'source':
                    qualifiers = feat.qualifiers
                    v['date'] = self.convert_gb_date(qualifiers['collection_date'][0])
                    v['country'] = qualifiers['country'][0]
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
