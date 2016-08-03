import os, re, time, datetime, csv, sys, json
import rethinkdb as r
from Bio import SeqIO
import argparse
from parse import parse
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='vdb', help="database to upload to")
parser.add_argument('--rethink_host', default=None, help="rethink host url")
parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
parser.add_argument('-v', '--virus', help="virus name")
parser.add_argument('--fname', help="input file name")
parser.add_argument('--ftype', default='fasta', help="input file format, default \"fasta\", other is \"genbank\", \"accession\" or \"tsv\"")
parser.add_argument('--path', default="data/", help="path to fasta file, default is \"data/\"")
parser.add_argument('--accessions', default=None, help="comma seperated list of accessions to be uploaded")
parser.add_argument('--email', default=None, help="email to access NCBI database via entrez to get virus information")
parser.add_argument('--overwrite', default=False, action="store_true",  help ="Overwrite fields that are not none")
parser.add_argument('--preview', default=False, action="store_true",  help ="If included, preview a virus document to be uploaded")
parser.add_argument('--replace', default=False, action="store_true",  help ="If included, delete all documents in table")

parser.add_argument('--source', default=None, help="source of fasta file")
parser.add_argument('--locus', default=None, help="gene or genomic region for sequences")
parser.add_argument('--host', default='human', help="host virus isolated from")
parser.add_argument('--country', default=None, help="country virus isolated from")
parser.add_argument('--authors', default=None, help="authors of source of sequences")
parser.add_argument('--private', default=False, action="store_true",  help ="sequences classified as not public")

class upload(parse):
    def __init__(self, database, virus, **kwargs):
        parse.__init__(self, **kwargs)
        self.virus = virus.lower()
        self.viruses_table = virus + "_viruses"
        self.sequences_table = virus + "_sequences"
        self.database = database.lower()
        self.uploadable_databases = ['vdb', 'test_vdb', 'test']
        if self.database not in self.uploadable_databases:
            raise Exception("Cant upload to this database: " + self.database, "add to list of databases allowed", self.uploadable_databases)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.viruses_table)
        self.rethink_io.check_table_exists(self.database, self.sequences_table)
        self.strains = {}

        self.virus_to_sequence_transfer_fields = []

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to VDB")
        viruses, sequences = self.parse(**kwargs)
        print('Formatting documents for upload')
        viruses = self.filter(viruses, **kwargs)
        sequences = self.filter(sequences, **kwargs)
        self.format(viruses, **kwargs)
        self.format(sequences, exclude_virus_methods=True, **kwargs)
        self.link_viruses_to_sequences(viruses, sequences)
        #self.transfer_fields(viruses, sequences, self.virus_to_sequence_transfer_fields)
        print("Upload Step")
        if not preview:
            print("Uploading viruses to " + self.database + "." + self.viruses_table)
            self.upload_documents(self.viruses_table, viruses, 'strain', **kwargs)
            print("Uploading sequences to " + self.database + "." + self.sequences_table)
            self.upload_documents(self.sequences_table, sequences, 'accession', **kwargs)
        else:
            print("Viruses:")
            print(json.dumps(viruses[0:2], indent=1))
            print("Sequences:")
            print(json.dumps(sequences[0:2], indent=1))
            print("Remove \"--preview\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def format(self, documents, exclude_virus_methods=False, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        self.define_regions()
        for doc in documents:
            if 'strain' in doc:
                doc['strain'] = self.fix_name(doc['strain'])
            self.format_date(doc)
            self.format_place(doc)
            self.format_region(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)
        if not exclude_virus_methods:
            print("Determining latitudes and longitudes")
            self.determine_latitude_longitude(documents, ['country'])

    def fix_name(self, name):
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '').replace('duck', '').replace('environment', '')
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
        date_fields = []
        for f in ['date', 'collection_date', 'submission_date']:
            if f in virus:
                date_fields.append(f)

        for field in date_fields:
            if virus[field] is not None:
                virus[field] = re.sub(r'_', r'-', virus[field])
                # ex. 2002 (Month and day unknown)
                if re.match(r'\d\d\d\d-(\d\d|XX)-(\d\d|XX)', virus[field]):
                    pass
                elif re.match(r'\d\d\d\d\s\(Month\sand\sday\sunknown\)', virus[field]):
                    virus[field] = virus[field][0:4] + "-XX-XX"
                # ex. 2009-06 (Day unknown)
                elif re.match(r'\d\d\d\d-\d\d\s\(Day\sunknown\)', virus[field]):
                    virus[field] = virus[field][0:7] + "-XX"
                elif re.match(r'\d\d\d\d-\d\d', virus[field]):
                    virus[field] = virus[field][0:7] + "-XX"
                elif re.match(r'\d\d\d\d', virus[field]):
                    virus[field] = virus[field][0:4] + "-XX-XX"
                else:
                    print("Couldn't reformat this date: " + virus[field] + ", setting to None")
                    virus[field] = None

    def camelcase_to_snakecase(self, name):
        '''
        convert camelcase format to snakecase format
        :param name:
        :return:
        '''
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

    def snakecase_to_camelcase(self, name):
        split_name = name.split('_')
        split_name = [x.title() for x in split_name]
        return "".join(split_name)

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
        if 'country' in virus:
            virus['region'] = '?'
            if virus['country'] is not None:
                test_country = self.camelcase_to_snakecase(virus['country'])
                if test_country is not None and test_country in self.country_to_region:
                    virus['region'] = self.country_to_region[self.camelcase_to_snakecase(virus['country'])]
                if virus['country'] != '?' and virus['region'] == '?':
                    virus['region'] = None
                    print("couldn't parse region for " + virus['strain'] + ", country: " + str(virus["country"]))

    def format_place(self, virus):
        '''
        Ensure snakecaase formatting after assigning these fields
        '''
        location_fields = ['country', 'division', 'location']
        for field in location_fields:
            if field in virus and virus[field] is not None:
                if "_" in virus[field]:  # French_Polynesia -> french_polynesia
                    virus[field] = "_".join(virus[field].split("_")).lower()
                virus[field] = self.camelcase_to_snakecase(virus[field])

    def filter(self, documents, **kwargs):
        return documents

    def determine_latitude_longitude(self, documents, location_fields=['location', 'division', 'country']):
        '''
        Assign latitude, longitude for as specific of a location as possible
        Location fields should be defined in order prioritized to determine latitude and longitude from
        Start by looking first for location, then division then country

        '''
        # get the latitude and longitudes that were already determined
        file = open("source-data/geo_lat_long.tsv", 'r')
        reader = csv.DictReader(filter(lambda row: row[0]!='#', file), delimiter='\t')		# list of dicts
        self.location_to_lat_long = {}
        for line in reader:
            self.location_to_lat_long[line['location'] + ":" + line['country_code']] = (float(line['latitude']), float(line['longitude']))
        file.close()

        # Mapping from country to ISO 3166-1 Alpha-2 code
        file = open("source-data/geo_ISO_code.tsv", 'r')
        reader = csv.DictReader(filter(lambda row: row[0]!='#', file), delimiter='\t')
        self.country_to_code = {}
        for line in reader:
            self.country_to_code[line['country']] =line['code']
        file.close()

        # file to write the new latitudes and longitudes that will be determined
        file = open("source-data/geo_lat_long.tsv", 'a')
        for doc in documents:
            if 'country' in doc:
                if doc['country'] in self.country_to_code:
                    country_code = self.country_to_code[doc['country']]
                    locations = [doc[field] for field in location_fields]
                    determined_location = False
                    for loc in locations:
                        if loc + ":" + country_code in self.location_to_lat_long:  # already determined location
                            doc['latitude'], doc['longitude'], doc['lat_long_location'] = self.location_to_lat_long[loc + ":" + country_code] + (loc,)
                            determined_location = True
                        else:
                            try:
                                result = self.get_latitude_longitude(country_code, loc)
                            except:
                                print(doc['strain'], locations, loc)
                                file.close()
                                raise Exception("Geopy query failed")
                            if result is not None:  # found location
                                if loc not in self.location_to_lat_long:
                                    file.write("\t".join([loc, country_code, str(result[0]), str(result[1])]) + "\n")
                                    self.location_to_lat_long[loc + ":" + country_code] = (result[0], result[1])
                                doc['latitude'], doc['longitude'], doc['lat_long_location'] = result + (loc,)
                                determined_location = True
                        if determined_location:
                            break
                    if not determined_location:
                        print("Couldn't determine latitude and longitude for ", doc['strain'], locations)
                        doc['latitude'], doc['longitude'], doc['lat_long_location'] = None, None, None
                else:
                    print("Couldn't find alpha-2 country code for ", doc['country'])
            else:
                print("Country not defined document, can't determine latitude and longitude", doc['strain'])

    def get_latitude_longitude(self, country_code, location):
        '''
        Use geopy package to determine latitude and longitude for location
        Bias results to country_code
        Return tuple of latitude, longitude if result found, otherwise returns None
        '''
        from geopy.geocoders import Nominatim
        geolocator = Nominatim(country_bias=country_code)
        result = geolocator.geocode(location.replace("_", " "))
        if result is not None:
            print("Found latitude and longitude for ", location)
            return (result.latitude, result.longitude)
        else:
            print("Couldn't find latitude and longitude for ", location)
            return result

    def link_viruses_to_sequences(self, viruses, sequences):
        '''
        Link the sequence information virus isolate information via the strain name
        '''
        strain_name_to_virus_doc = {virus['strain']: virus for virus in viruses}
        for sequence_doc in sequences:
            if sequence_doc['strain'] in strain_name_to_virus_doc:  # determine if sequence has a corresponding virus to link to
                virus_doc = strain_name_to_virus_doc[sequence_doc['strain']]
                virus_doc['sequences'].append(sequence_doc['accession'])
                virus_doc['number_sequences'] += 1

    def transfer_fields(self, docs_from, docs_to, fields=[]):
        '''
        transfer fields between corresponding dictionaries based on matching strain field
        '''
        if len(fields) > 0:
            strain_name_to_docs_from = {doc['strain']: doc for doc in docs_from}
            for field in fields:
                for doc_to in docs_to:
                    strain = doc_to['strain']
                    doc_to[field] = strain_name_to_docs_from[strain][field]
        # delete fields after adding to every doc_to
        for doc_from in docs_from:
            for field in fields:
                del doc_from[field]

    def upload_documents(self, table, documents, index, replace, **kwargs):
        if replace:
            print("Deleting documents in database:" + self.database + "." + table)
            r.table(table).delete().run()
        db_documents = list(r.db(self.database).table(table).run())
        db_relaxed_keys = self.relaxed_keys(db_documents)
        db_key_to_documents = {dd[index]: dd for dd in db_documents}
        update_documents = {}
        upload_documents = {}
        print("Classifying documents to upload or update")
        # Determine whether documents need to be updated or uploaded
        for doc in documents:
            # match to relaxed strain name if available
            db_key = doc[index]
            if self.relax_name(doc[index]) in db_relaxed_keys:
                db_key = db_relaxed_keys[self.relax_name(doc[index])]

            if db_key in db_key_to_documents.keys():  # add to update documents
                update_documents[db_key] = doc
            elif db_key in upload_documents.keys():  # document already in list to be uploaded, check for updates
                self.update_document_meta(upload_documents[db_key], doc, output=False, **kwargs)
            else:  # add to upload documents
                upload_documents[doc[index]] = doc
        print("Inserting ", len(upload_documents), "documents")
        self.upload_to_rethinkdb(self.database, table, upload_documents.values(), 'error')
        self.check_for_updates(table, update_documents, db_key_to_documents, **kwargs)

    def upload_to_rethinkdb(self, database, table, documents, conflict_resolution):
        optimal_upload = 200
        if len(documents) > optimal_upload:
            list_documents = [documents[x:x+optimal_upload] for x in range(0, len(documents), optimal_upload)]
        else:
            list_documents = [documents]
        print("Uploading to rethinkdb in " + str(len(list_documents)) + " batches")
        for list_docs in list_documents:
            try:
                r.table(table).insert(list_docs, conflict=conflict_resolution).run()
            except:
                raise Exception("Couldn't insert new documents into database", database + "." + table)

    def check_for_updates(self, table, update_documents, db_key_to_documents, **kwargs):
        print("Checking for updates to ", len(update_documents), "documents")
        # determine which documents need to be updated and update them
        updated = [db_key_to_documents[db_key] for db_key, doc in update_documents.items() if self.update_document_meta(db_key_to_documents[db_key], doc, **kwargs)]
        # then insert the updated documents
        if len(updated) > 0:
            print("Found updates to ", len(updated), "documents")
            self.upload_to_rethinkdb(self.database, table, updated, 'replace')
        else:
            print("No documents need to be updated in ", self.database + "." + table)

    def update_document_meta(self, db_doc, doc, overwrite, output=True, **kwargs):
        '''
        update overwritable fields at the base level of the document
        Updates the db_doc to fields of doc based on rules below
        '''
        updated = False
        for field in doc:
            if field == 'timestamp':
                pass
            elif field == 'sequences':
                # add new accessions to sequences
                for accession in doc[field]:
                    if accession not in db_doc[field]:
                        db_doc[field].append(accession)
                        db_doc['number_sequences'] += 1
            else:
                # update if field not present in db_doc
                if field not in db_doc:
                    if doc[field] is not None:
                        print("Creating field ", field, " assigned to ", doc[field])
                    db_doc[field] = doc[field]
                    updated = True
                #update db_doc information if doc info is different, or if overwrite is false update db_doc information only if doc info is different and not null
                elif doc[field] is not None and (overwrite and db_doc[field] != doc[field]) or (not overwrite and db_doc[field] is None and db_doc[field] != doc[field]):
                    if output:
                        print("Updating field " + str(field) + ", from \"" + str(db_doc[field]) + "\" to \"" + str(doc[field])) + "\""
                    db_doc[field] = doc[field]
                    updated = True
        return updated

    def relaxed_keys(self, documents):
        '''
        Create dictionary from relaxed vdb strain names to actual vdb strain names.
        '''
        strains = {}
        for doc in documents:
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

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0:'accession', 1:'strain', 2:'date', 4:'country', 5:'division', 6:'location'}
    setattr(args, 'fasta_fields', fasta_fields)
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
