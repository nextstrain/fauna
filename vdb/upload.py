import os, re, time, datetime, csv, sys, json
import rethinkdb as r
from Bio import SeqIO
import argparse
from parse import parse
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to upload to")
parser.add_argument('-v', '--virus', help="virus table to interact with")
parser.add_argument('--fname', help="input file name")
parser.add_argument('--ftype', default='fasta', help="input file format, default \"fasta\", other is \"genbank\" or \"accession\"")
parser.add_argument('--accessions', default=None, help="comma seperated list of accessions to be uploaded")
parser.add_argument('--source', default=None, help="source of fasta file")
parser.add_argument('--locus', default=None, help="gene or genomic region for sequences")
parser.add_argument('--host', default='human', help="host virus isolated from")
parser.add_argument('--authors', default=None, help="authors of source of sequences")
parser.add_argument('--private', default=False, action="store_true",  help ="sequences classified as not public")
parser.add_argument('--path', default="data/", help="path to fasta file, default is \"data/\"")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('--preview', default=False, action="store_true",  help ="If included, preview a virus document to be uploaded")
parser.add_argument('--overwrite', default=False, action="store_true",  help ="Overwrite fields that are not none")
parser.add_argument('--exclusive', default=True, action="store_false",  help ="download all docs in db to check before upload")
parser.add_argument('--email', default=None, help="email to access NCBI database via entrez to get virus information")
parser.add_argument('--auto_upload', default=False, action="store_true", help="search genbank for recent sequences")

class upload(parse):
    def __init__(self, database, virus, **kwargs):
        parse.__init__(self, **kwargs)
        self.virus = virus.lower()
        self.database = database.lower()
        self.uploadable_databases = ['vdb', 'test_vdb', 'test']
        if self.database not in self.uploadable_databases:
            raise Exception("Cant upload to this database: " + self.database, "add to list of databases allowed", self.uploadable_databases)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.virus)

        # fields that are needed to upload
        self.index_field = ['strain']
        self.sequence_upload_fields = ['locus', 'sequence']
        self.sequence_optional_fields = ['accession']  # ex. if from virological.org or not in a database
        self.citation_upload_fields = ['source']
        self.citation_optional_fields = ['authors', 'title', 'url']
        self.grouping_upload_fields = []  # need to define for each virus specific upload script if applicable
        self.grouping_optional_fields = []
        self.virus_upload_fields = ['strain', 'date', 'country', 'virus', 'date_modified', 'public', 'region', 'host'] # sequences, citations, added after filtering
        self.virus_optional_fields = ['division', 'location']
        self.upload_fields = self.virus_upload_fields+self.sequence_upload_fields+self.citation_upload_fields+self.grouping_upload_fields
        self.optional_fields = self.virus_optional_fields+self.sequence_optional_fields+self.citation_optional_fields+self.grouping_optional_fields
        self.overwritable_virus_fields = ['date', 'country', 'division', 'location', 'virus', 'public']
        self.strains = {}

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to VDB")
        self.parse(**kwargs)
        self.format()
        self.filter()
        if not preview:
            self.upload_documents(**kwargs)
        else:
            try:
                print(json.dumps(self.viruses[0], indent=1))
            except:
                print(json.dumps(self.viruses, indent=1))
            print("Include \"--upload\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def format(self):
        '''
        format virus information in preparation to upload to database table
        '''
        print('Formatting for upload')
        self.define_regions()
        for virus in self.viruses:
            self.format_date(virus)
            self.format_place(virus)
            self.format_region(virus)
            self.rethink_io.delete_extra_fields(virus, self.upload_fields+self.optional_fields, self.index_field)

    def fix_name(self, name):
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '')
        try:
            tmp_name = 'V' + str(int(tmp_name))
        except:
            pass
        return tmp_name

    def format_date(self, virus):
        '''
        Format viruses date attribute: collection date in YYYY-MM-DD format, for example, 2016-02-28
        Input date could be YYYY_MM_DD, reformat to YYYY-MM-DD
        '''
        # ex. 2002_04_25 to 2002-04-25
        virus['date'] = re.sub(r'_', r'-', virus['date'])
        # ex. 2002 (Month and day unknown)
        if re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', virus['date']):
            pass
        elif re.match(r'\d\d\d\d\s\(Month\sand\sday\sunknown\)', virus['date']):
            virus['date'] = virus['date'][0:4] + "-XX-XX"
        # ex. 2009-06 (Day unknown)
        elif re.match(r'\d\d\d\d-\d\d\s\(Day\sunknown\)', virus['date']):
            virus['date'] = virus['date'][0:7] + "-XX"
        elif re.match(r'\d\d\d\d-\d\d', virus['date']):
            virus['date'] = virus['date'][0:7] + "-XX"
        elif re.match(r'\d\d\d\d', virus['date']):
            virus['date'] = virus['date'][0:4] + "-XX-XX"
        else:
            print("Couldn't reformat this date: " + virus['date'])

    def camelcase_to_snakecase(self, name):
        '''
        convert camelcase format to snakecase format
        :param name:
        :return:
        '''
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

    def define_regions(self):
        '''
        open country to region dictionary
        '''
        try:
            reader = csv.DictReader(open("source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
        except:
            raise Exception("Couldn't find geo regions file")
        self.country_to_region = {}
        for line in reader:
            self.country_to_region[line['country']] = line['region']

    def format_region(self, virus):
        '''
        Label viruses with region based on country, if prune then filter out viruses without region
        '''
        virus['region'] = '?'
        test_country = virus['country']
        if test_country is not None and test_country in self.country_to_region:
            virus['region'] = self.country_to_region[virus['country']].lower()
        if virus['country'] != '?' and virus['region'] == '?':
            virus['region'] = None
            print("couldn't parse region for " + virus['strain'] + ", country: " + str(virus["country"]))

    def format_place(self, virus):
        '''
        Ensure lowercase_with_underscores formatting after assigning these fields
        '''
        location_fields = ['region', 'country', 'division', 'location']
        for field in location_fields:
            if field in virus and virus[field] is not None:
                virus[field] = virus[field].lower().replace(' ', '_')

    def filter(self):
        '''
        Filter out viruses without correct format, check  optional and upload attributes
        '''
        print(str(len(self.viruses)) + " viruses before filtering")
        self.rethink_io.check_optional_attributes(self.viruses, self.optional_fields)
        self.viruses = filter(lambda v: re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', v['date']), self.viruses)
        self.viruses = filter(lambda v: isinstance(v['public'], (bool)), self.viruses)
        self.viruses = filter(lambda v: v['region'] is not None, self.viruses)
        self.viruses = filter(lambda v: self.rethink_io.check_required_attributes(v, self.upload_fields, self.index_field), self.viruses)
        print(str(len(self.viruses)) + " viruses after filtering")
        self.format_schema()

    def format_schema(self):
        '''
        move sequence information into nested 'sequences' field
        '''
        for virus in self.viruses:
            virus['sequences'] = [{}]
            virus['citations'] = [{}]
            for field in self.sequence_upload_fields + self.sequence_optional_fields:
                if field in virus.keys():
                    virus['sequences'][0][field] = virus[field]
                    del virus[field]
            for field in self.citation_upload_fields + self.citation_optional_fields:
                if field in virus.keys():
                    virus['citations'][0][field] = virus[field]
                    del virus[field]

    def upload_documents(self, exclusive, **kwargs):
        '''
        Insert viruses into collection
        '''
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        db_relaxed_strains = self.relaxed_strains()
        # Faster way to upload documents, downloads all database documents locally and looks for precense of strain in database
        if exclusive:
            db_viruses = list(r.db(self.database).table(self.virus).run())
            db_strain_to_viruses = {db_v['strain']: db_v for db_v in db_viruses}
            update_viruses = {}
            upload_viruses = {}
            for virus in self.viruses:
                # determine the corresponding database strain name based on relaxed db and virus strain
                db_strain = virus['strain']
                if self.relax_name(virus['strain']) in db_relaxed_strains:
                    db_strain = db_relaxed_strains[self.relax_name(virus['strain'])]
                if db_strain in db_strain_to_viruses.keys():  # virus already in database
                    update_viruses[db_strain] = virus
                elif db_strain in upload_viruses.keys():  # virus already to be uploaded, need to check for updates to sequence information
                    upload_v = upload_viruses[db_strain]
                    self.update_document_sequence(upload_v, virus, **kwargs)  # add new sequeunce information to upload_v
                else:  # new virus that needs to be uploaded
                    upload_viruses[virus['strain']] = virus
            print("Inserting ", len(upload_viruses), "viruses into database", self.virus)
            try:
                r.table(self.virus).insert(upload_viruses.values()).run()
            except:
                raise Exception("Couldn't insert new viruses into database")
            print("Checking for updates to ", len(update_viruses), "viruses in database", self.virus)
            updated = []
            for db_strain, v in update_viruses.items():  # determine if virus has new information
                document = db_strain_to_viruses[db_strain]
                updated_sequence = self.update_document_sequence(document, v, **kwargs)
                updated_meta = self.update_document_meta(document, v, self.overwritable_virus_fields, **kwargs)
                if updated_sequence or updated_meta:
                    document['date_modified'] = v['date_modified']
                    updated.append(document)
            try:
                r.table(self.virus).insert(updated, conflict="replace").run()
            except:
                raise Exception("Couldn't update viruses already in database")
        # Slower way to upload but less chance database is changed while updating documents, asks database for each document separately
        else:
            print("Uploading " + str(len(self.viruses)) + " viruses to the table")
            self.relaxed_strains()
            for virus in self.viruses:
                # Retrieve virus from table to see if it already exists, try relaxed comparison first
                relaxed_name = virus['strain']
                if self.relax_name(virus['strain']) in db_relaxed_strains:
                    relaxed_name = db_relaxed_strains[self.relax_name(virus['strain'])]
                try:
                    document = r.table(self.virus).get(relaxed_name).run()
                except:
                    print(virus)
                    raise Exception("Couldn't retrieve this virus")
                # Virus doesn't exist in table yet so add it
                if document is None:
                    #print("Inserting " + virus['strain'] + " into database")
                    try:
                        r.table(self.virus).insert(virus).run()
                    except:
                        print(virus)
                        raise Exception("Couldn't insert this virus")
                # Virus exists in table so just add sequence information and update meta data if needed
                else:
                    updated_sequence = self.update_document_sequence(document, virus, **kwargs)
                    updated_meta = self. update_document_meta(relaxed_name, document, virus, self.overwritable_virus_fields, **kwargs)
                    if updated_sequence or updated_meta:
                        document['date_modified'] = virus['date_modified']
                        r.table(self.virus).insert(document, conflict="replace").run()

    def relaxed_strains(self):
        '''
        Create dictionary from relaxed vdb strain names to actual vdb strain names.
        '''
        strains = {}
        cursor = list(r.db(self.database).table(self.virus).run())
        for doc in cursor:
            strains[self.relax_name(doc['strain'])] = doc['strain']
        return strains

    def relax_name(self, name):
        '''
        Return the relaxed strain name to compare with
        '''
        name = re.sub(r"-", '', name)
        name = re.sub(r"_", '', name)
        name = re.sub(r"/", '', name)
        return name

    def update_document_meta(self, document, v, overwritable_fields, overwrite, **kwargs):
        '''
        update overwritable fields at the base level of the document
        '''
        updated = False
        for field in overwritable_fields:
            # update if field not present in document
            if field not in document:
                if field in v:
                    print("Creating virus field ", field, " assigned to ", v[field])
                    document[field] = v[field]
                    updated = True
            #update doc_virus information if virus info is different, or if overwrite is false update doc_virus information only if virus info is different and not null
            elif (overwrite and document[field] != v[field]) or (not overwrite and document[field] is None and document[field] != v[field]):
                if field in v:
                    print("Updating virus field " + str(field) + ", from \"" + str(document[field]) + "\" to \"" + v[field]) + "\""
                    document[field] = v[field]
                    updated = True
        return updated

    def update_document_sequence(self, document, v, **kwargs):
        '''
        Update sequence fields if matching accession or sequence
        Append sequence to sequence list if no matching accession or sequence
        '''
        doc_seqs = document['sequences']
        virus_seq = v['sequences'][0]
        virus_citation = v['citations'][0]
        # try comparing accession's first
        if virus_seq['accession'] != None:
            if all(virus_seq['accession'] != seq_info['accession'] for seq_info in doc_seqs):
                updated = self.append_new_sequence(document, virus_seq, virus_citation)
            else:
                updated = self.update_sequence_citation_field(document, v, 'accession', self.sequence_upload_fields+self.sequence_optional_fields, self.citation_optional_fields, **kwargs)
        else:  # if it doesn't have an accession, see if sequence already in document
            if all(virus_seq['sequence'] != seq_info['sequence'] for seq_info in doc_seqs):
                updated = self.append_new_sequence(document, virus_seq, virus_citation)
            else:
                updated = self.update_sequence_citation_field(document, v, 'sequence', self.sequence_upload_fields+self.sequence_optional_fields, self.citation_optional_fields, **kwargs)
        if len(document['sequences']) != len(document['citations']):
            print("Warning: length of list of sequences and citations does not match for " + v['strain'])
        return updated

    def append_new_sequence(self,document, virus_seq, virus_citation):
        '''
        New sequence information, so append to document sequences and citations field
        '''
        document['sequences'].append(virus_seq)
        document['citations'].append(virus_citation)
        return True

    def update_sequence_citation_field(self, document, virus_doc, check_field, sequence_fields, citation_fields, **kwargs):
        '''
        Based on the check_field (accession or sequence), look for updated sequence and citation information
        '''
        updated_sequence = False
        updated_citation = False
        strain = virus_doc['strain']
        doc_seqs = document['sequences']
        virus_seq = virus_doc['sequences'][0]
        index = -1
        for doc_sequence_info in doc_seqs:
            index += 1
            if doc_sequence_info[check_field] == virus_seq[check_field]:  # find the identical sequence info based on check field
                updated_sequence = self.update_nested_field(strain, sequence_fields, doc_sequence_info, virus_seq, **kwargs)
                updated_citation = self.update_nested_field(strain, citation_fields, document['citations'][index], virus_doc['citations'][0], **kwargs)
        updated = False
        if updated_sequence:
            document['sequences'] = doc_seqs
            updated = True
        if updated_citation:
            document['citations'][index] = document['citations'][index]
            updated = True
        return updated

    def update_nested_field(self, strain, fields, doc_info, virus_info, overwrite, **kwargs):
        '''
        check for updates to nested sequences and citations fields.
        '''
        updated_sequence = False
        for field in fields:
            if field not in doc_info:
                if field in virus_info:
                    print("Creating field ", field, " assigned to \"", virus_info[field] + "\" for strain: " + strain)
                    doc_info[field] = virus_info[field]
                    updated_sequence = True
            elif (overwrite and doc_info[field] != virus_info[field]) or (not overwrite and doc_info[field] is None and doc_info[field] != virus_info[field]):
                if field in virus_info:
                    print("Updating field " + str(field) + ", from \"" + str(doc_info[field]) + "\" to \"" + str(virus_info[field]) + "\" for strain: " + strain)
                    doc_info[field] = virus_info[field]
                    updated_sequence = True
        return updated_sequence

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 1:'strain', 2:'date', 4:'country', 5:'division', 6:'location'}
    setattr(args, 'fasta_fields', fasta_fields)
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
