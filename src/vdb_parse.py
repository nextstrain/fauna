import os, re, datetime
from Bio import SeqIO

class vdb_parse(object):
    def __init__(self, fasta_fields, **kwargs):
        self.fasta_fields = fasta_fields

    def parse(self):
        if self.ftype == 'genbank':
            self.viruses = self.parse_gb(self.path + self.fasta_fname)
        elif self.ftype == 'fasta':
            self.viruses = self.parse_fasta(self.path + self.fasta_fname)
        return self.viruses

    def parse_fasta(self, fasta):
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

    def parse_gb(self, gb):
        '''
        Parse genbank file
        :return: list of documents(dictionaries of attributes) to upload
        '''
        viruses = []
        try:
            handle = open(gb, 'r')
        except IOError:
            print(gb, "not found")
        else:
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
                if reference.title is not None:
                    v['title'] = reference.title
                v['url'] = "http://www.ncbi.nlm.nih.gov/nuccore/" + v['accession']
                record_features = record.features
                for feat in record_features:
                    if feat.type == 'source':
                        qualifiers = feat.qualifiers
                        v['date'] = self.convert_gb_date(qualifiers['collection_date'][0])
                        v['country'] = qualifiers['country'][0]
                        v['strain'] = qualifiers['isolate'][0]
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
        Credit to Gytis Dudas
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
