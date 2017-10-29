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
parser.add_argument('--fasta_header_fix', default=None, help="faster header fix file. Default: None.")
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
parser.add_argument('--title', default=None, help="title of sequence release")
parser.add_argument('--url', default=None, help="url of source of sequences")
parser.add_argument('--public', default=True, dest='public', action="store_true",  help ="sequences classified as public")
parser.add_argument('--private', default=False, dest='public', action="store_false",  help ="sequences classified as private")

class upload(parse):
    def __init__(self, database, virus, **kwargs):
        parse.__init__(self, **kwargs)
        self.virus = virus.lower()
        self.viruses_table = virus + "_viruses"
        self.sequences_table = virus + "_sequences"
        self.database = database.lower()
        self.uploadable_databases = ['vdb', 'test_vdb', 'test']
        self.strains = {}
        self.strain_fix_fname = None
        self.location_fix_fname = None
        self.fix_location = None
        self.date_fix_fname = None
        self.fix_date = None
        self.virus_to_sequence_transfer_fields = []

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        self.connect(**kwargs)
        print("Uploading Viruses to VDB")
        viruses, sequences = self.parse(**kwargs)
        print('Formatting documents for upload')
        self.format_viruses(viruses, **kwargs)
        self.format_sequences(sequences, **kwargs)
        print("")
        print("Filtering out viruses")
        viruses = self.filter(viruses, 'strain', **kwargs)
        print("Filtering out sequences")
        sequences = self.filter(sequences, 'accession', **kwargs)
        print("")
        #self.match_duplicate_strains(viruses, sequences, **kwargs)
        #self.match_database_duplicate_strains(viruses, sequences, **kwargs)
        #self.match_duplicate_accessions(sequences, **kwargs)
        #self.match_database_duplicate_accessions(sequences, **kwargs)
        self.link_viruses_to_sequences(viruses, sequences)
        #self.transfer_fields(viruses, sequences, self.virus_to_sequence_transfer_fields)
        print("")
        print("Upload Step")
        if not preview:
            print("Uploading viruses to " + self.database + "." + self.viruses_table)
            self.upload_documents(self.viruses_table, viruses, index='strain', **kwargs)
            print("Uploading sequences to " + self.database + "." + self.sequences_table)
            self.upload_documents(self.sequences_table, sequences, index='accession', **kwargs)
        else:
            print("Viruses:")
            print(json.dumps(viruses[0], indent=1))
            print("Sequences:")
            print(json.dumps(sequences[0], indent=1))
            print("Remove \"--preview\" to upload documents")
            print("Printed preview of viruses to be uploaded to make sure fields make sense")

    def connect(self, **kwargs):
        if self.database not in self.uploadable_databases:
            raise Exception("Can't upload to this database: " + self.database, "add to list of databases allowed", self.uploadable_databases)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.viruses_table)
        self.rethink_io.check_table_exists(self.database, self.sequences_table)

    def format_viruses(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        if self.strain_fix_fname is not None:
            self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        if self.location_fix_fname is not None:
            self.fix_location = self.define_location_fixes(self.location_fix_fname)
        if self.date_fix_fname is not None:
            self.fix_date = self.define_date_fixes(self.date_fix_fname)
        self.define_regions("source-data/geo_regions.tsv")
        self.define_countries("source-data/geo_synonyms.tsv")
        self.define_latitude_longitude("source-data/geo_lat_long.tsv", "source-data/geo_ISO_code.tsv")
        for doc in documents:
            if 'strain' not in doc:
                doc['strain'] = "unnamed"
            doc['strain'], doc['original_strain'] = self.fix_name(doc['strain'])
            if self.fix_location is not None:
                if doc['strain'] in self.fix_location:
                    doc['location'] = self.fix_location[doc['strain']]
            if self.fix_date is not None:
                if doc['strain'] in self.fix_date:
                    doc['collection_date'] = self.fix_date[doc['strain']]
            self.format_date(doc)
            self.format_place(doc)
            self.format_region(doc)
            self.determine_latitude_longitude(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)

    def format_sequences(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        if self.strain_fix_fname is not None:
            self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        for doc in documents:
            if 'strain' not in doc:
                doc['strain'] = "unnamed"
            doc['strain'], doc['original_strain'] = self.fix_name(doc['strain'])
            self.format_date(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)

    def define_strain_fixes(self, fname):
        '''
        Open strain name fixing files and define corresponding dictionaries
        '''
        reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
        fix_whole_name = {}
        for line in reader:
            fix_whole_name[line['label'].decode('unicode-escape')] = line['fix']
        return fix_whole_name

    def define_location_fixes(self, fname):
        '''
        Open location fix file and define corresponding dictionaries
        '''
        reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
        fix_location = {}
        for line in reader:
            fix_location[line['label'].decode('unicode-escape')] = line['fix']
        return fix_location

    def define_date_fixes(self, fname):
        '''
        Open date fix file and define corresponding dictionaries
        '''
        reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
        fix_date = {}
        for line in reader:
            fix_date[line['label'].decode('unicode-escape')] = line['fix']
        return fix_date

    def replace_strain_name(self, original_name, fixes={}):
        '''
        return the new strain name that will replace the original
        '''
        if original_name in fixes:
            return fixes[original_name]
        else:
            return original_name

    def fix_name(self, name):
        original_name = name
        name = original_name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '').replace('duck', '').replace('environment', '')
        try:
            name = 'V' + str(int(name))
        except:
            pass
        return name, original_name

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
            if virus[field] is not None and virus[field].strip() != '':
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
            else:
                virus[field] = None

    def camelcase_to_snakecase(self, name):
        '''
        convert camelcase format to snakecase format
        :param name:
        :return:
        '''
        if name is not None:
            s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
            return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower().replace(" ", "")

    def snakecase_to_camelcase(self, name):
        if name is not None:
            split_name = name.split('_')
            split_name = [x.title() for x in split_name]
            return "".join(split_name)

    def define_countries(self, fname):
        '''
        open synonym to country dictionary
        Location is to the level of country of administrative division when available
        '''
        reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')		# list of dicts
        self.label_to_location = {}
        self.label_to_division = {}
        self.label_to_country = {}

        for line in reader:
            self.label_to_location[line['label'].decode('unicode-escape').lower()] = self.camelcase_to_snakecase(line['location'])
            self.label_to_division[line['label'].decode('unicode-escape').lower()] = self.camelcase_to_snakecase(line['division'])
            self.label_to_country[line['label'].decode('unicode-escape').lower()] = self.camelcase_to_snakecase(line['country'])

    def define_regions(self, fname):
        '''
        open country to region dictionary
        '''
        try:
            reader = csv.DictReader(open(fname), delimiter='\t')		# list of dicts
        except:
            raise Exception("Couldn't find geo regions file")
        self.country_to_region = {}
        for line in reader:
            self.country_to_region[line['country']] = line['region']

    def format_place(self, doc, determine_location=True):
        '''
        Try to determine location information from location fields
        Ensure snakecase formatting after assigning location fields
        '''
        location_fields = ['location', 'division', 'country']
        for field in location_fields:
            if determine_location:
                if field in doc and doc[field] is not None:
                    result = self.determine_location(self.snakecase_to_camelcase(doc[field]))
                    if result is not None:
                        doc['location'], doc['division'], doc['country'] = result
                        break
                    else:
                        print("couldn't parse location for ", doc['strain'], self.snakecase_to_camelcase(doc[field]))
                        if "_" in doc[field]:  # French_Polynesia -> french_polynesia
                            doc[field] = "_".join(doc[field].split("_")).lower()
                        doc[field] = self.camelcase_to_snakecase(doc[field])
            else:
                if field in doc and doc[field] is not None:
                    if "_" in doc[field]:  # French_Polynesia -> french_polynesia
                        doc[field] = "_".join(doc[field].split("_")).lower()
                    doc[field] = self.camelcase_to_snakecase(doc[field])

    def determine_location(self, name):
        '''
        Try to determine country, division and location information from name
        Return tuple of country, division, location if found, otherwise return None
        '''
        try:
            if name in self.label_to_country:
                return (self.label_to_location[name], self.label_to_division[name], self.label_to_country[name])
            else:
                label = re.match(r'^([^/]+)', name).group(1).lower()						# check first for whole geo match
            if label in self.label_to_country:
                return ( self.label_to_location[label], self.label_to_division[label], self.label_to_country[label])
            else:
                label = re.match(r'^([^\-^\/]+)', name).group(1).lower()			# check for partial geo match A/CHIBA-C/61/2014
            if label in self.label_to_country:
                return ( self.label_to_location[label], self.label_to_division[label], self.label_to_country[label])
            else:
                label = re.match(r'^([A-Z][a-z]+)[A-Z0-9]', name).group(1).lower()			# check for partial geo match
            if label in self.label_to_country:
                return ( self.label_to_location[label], self.label_to_division[label], self.label_to_country[label])
            else:
                return None
        except:
            return None

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

    def filter(self, documents, index, **kwargs):
        '''
        filter out certain documents
        '''
        print(str(len(documents)) + " documents before filtering")
        documents = filter(lambda doc: index in doc, documents)
        print(str(len(documents)) + " documents after filtering")
        return documents

    def define_latitude_longitude(self, lat_long_fname, code_fname):
        # get the latitude and longitudes that were already determined
        file = open(lat_long_fname, 'r')
        reader = csv.DictReader(filter(lambda row: row[0]!='#', file), delimiter='\t')		# list of dicts
        self.location_to_lat_long = {}
        for line in reader:
            try:
                self.location_to_lat_long[line['location'] + ":" + line['country_code']] = (float(line['latitude']), float(line['longitude']))
            except:
                print("Line failed ", line)
                raise Exception("Failed to read ", file, "please check the line that failed")
        file.close()

        # Mapping from country to ISO 3166-1 Alpha-2 code
        file = open(code_fname, 'r')
        reader = csv.DictReader(filter(lambda row: row[0]!='#', file), delimiter='\t')
        self.country_to_code = {}
        for line in reader:
            self.country_to_code[line['country']] =line['code']
        file.close()

        # file to write the new latitudes and longitudes that will be determined
        self.new_lat_long_file = open(lat_long_fname, 'a')

    def determine_latitude_longitude(self, doc, location_fields=['location', 'division', 'country']):
        '''
        Assign latitude, longitude for as specific of a location as possible
        Location fields should be defined in order prioritized to determine latitude and longitude from
        Start by looking first for location, then division then country

        '''
        if 'country' in doc and doc['country'] is not None:
            if doc['country'] in self.country_to_code:
                country_code = self.country_to_code[doc['country']]
                locations = [doc[field] for field in location_fields if field in doc and doc[field] is not None]
                for loc in locations:
                    if loc + ":" + country_code in self.location_to_lat_long:  # already determined location
                        doc['latitude'], doc['longitude'], doc['lat_long_location'] = self.location_to_lat_long[loc + ":" + country_code] + (loc,)
                        return doc['latitude'], doc['longitude'], doc['lat_long_location']
                    else:
                        try:
                            result = self.get_latitude_longitude(country_code, loc)
                        except:
                            print(doc['strain'], locations, loc)
                            self.new_lat_long_file.close()
                            raise Exception("Geopy query failed")
                        if result is not None:  # found location
                            if loc not in self.location_to_lat_long:
                                print("Found latitude and longitude for ", loc, doc['country'])
                                self.new_lat_long_file.write("\t".join([loc, country_code, str(result[0]), str(result[1])]) + "\n")
                                self.location_to_lat_long[loc + ":" + country_code] = (result[0], result[1])
                            doc['latitude'], doc['longitude'], doc['lat_long_location'] = result + (loc,)
                            return doc['latitude'], doc['longitude'], doc['lat_long_location']
                        else:
                            print("Couldn't find latitude and longitude for ", loc, doc['country'])
                print("Couldn't determine latitude and longitude for ", doc['strain'], locations)
                doc['latitude'], doc['longitude'], doc['lat_long_location'] = None, None, None
                return None, None, None
            else:
                print("Couldn't find alpha-2 country code for ", doc['country'], doc['strain'])
        else:
            #print("Country not defined for this document, can't determine latitude and longitude", doc['strain'])
            pass

    def get_latitude_longitude(self, country_code, location):
        '''
        Use geopy package to determine latitude and longitude for location
        Bias results to country_code
        Return tuple of latitude, longitude if result found, otherwise returns None
        '''
        from geopy.geocoders import Nominatim
        geolocator = Nominatim(country_bias=country_code)
        #from geopy.geocoders import GeoNames
        #geolocator = GeoNames(country_bias=country_code, username='')
        result = geolocator.geocode(location.replace("_", " "))
        if result is not None:
            return (result.latitude, result.longitude)
        else:
            return result

    def match_duplicate_strains(self, viruses, sequences, **kwargs):
        '''
        Compare strain names to other documents in the current upload to adjust their strain names to be the same
        Adjust the sequence strain name as well
        '''
        print("Adjusting strain names to match identical strains in documents to be uploaded")
        indexes = [doc['strain'] for doc in viruses]
        relaxed_strains = self.relaxed_keys(indexes, self.relax_name)
        adjusted_virus_strains = 0
        adjusted_sequence_strains = 0
        for doc in viruses:
            doc['strain'], adjusted_virus_strains = self.adjust_name(doc['strain'], relaxed_strains, adjusted_virus_strains)
        for doc in sequences:
            doc['strain'], adjusted_sequence_strains = self.adjust_name(doc['strain'], relaxed_strains, adjusted_sequence_strains)
        if adjusted_virus_strains > 0 or adjusted_sequence_strains > 0:
            print(str(adjusted_virus_strains) + " out of " + str(len(viruses)) + " virus strains were adjusted to match a virus in current upload documents")
            print(str(adjusted_sequence_strains) + " out of " + str(len(sequences)) + " sequence strains were adjusted to match a virus in current upload documents")

    def match_duplicate_accessions(self, sequences, **kwargs):
        '''
        Compare accessions to other documents in the current upload to adjust their accessions to be the same
        '''
        print("Adjusting accessions to match identical sequences in documents to be uploaded")
        indexes = [doc['accession'] for doc in sequences]
        relaxed_accessions = self.relaxed_keys(indexes, self.relax_name)
        adjusted_sequence_accessions = 0
        for doc in sequences:
            doc['accession'], adjusted_sequence_strains = self.adjust_name(doc['accession'], relaxed_accessions, adjusted_sequence_accessions)
        if adjusted_sequence_accessions > 0:
            print(str(adjusted_sequence_accessions) + " out of " + str(len(sequences)) + " sequence accessions were adjusted to match a sequence in current upload documents")

    def match_database_duplicate_strains(self, viruses, sequences, virus, database='vdb', **kwargs):
        '''
        Compare measurement strain names to database strain names to matching to existing viruses sequences
        '''
        table = virus + "_viruses"
        print("Using " + database + "." + table + " to adjust strain names to match strains already in " + database + "." + table)
        vdb_strains = self.relaxed_keys(set(list(r.db(database).table(table).get_field('strain').run())), self.relax_name)
        adjusted_virus_strains = 0
        adjusted_sequence_strains = 0
        for doc in viruses:
            doc['strain'], adjusted_virus_strains = self.adjust_name(doc['strain'], vdb_strains, adjusted_virus_strains)
        for doc in sequences:
            doc['strain'], adjusted_sequence_strains = self.adjust_name(doc['strain'], vdb_strains, adjusted_sequence_strains)
        if adjusted_virus_strains > 0 or adjusted_sequence_strains > 0:
            print(str(adjusted_virus_strains) + " out of " + str(len(viruses)) + " virus strains were adjusted to match a virus in " + database + "." + table)
            print(str(adjusted_sequence_strains) + " out of " + str(len(sequences)) + " sequence strains were adjusted to match a virus in " + database + "." + table)

    def match_database_duplicate_accessions(self, sequences, virus, database='vdb', **kwargs):
        '''
        Compare measurement strain names to database strain names to matching to existing viruses sequences
        '''
        table = virus + "_sequences"
        print("Using " + database + "." + table + " to adjust accessions to match sequences already in " + database + "." + table)
        vdb_strains = self.relaxed_keys(set(list(r.db(database).table(table).get_field('accession').run())), self.relax_name)
        adjusted_sequence_accessions = 0
        for doc in sequences:
            doc['accession'], adjusted_sequence_accessions = self.adjust_name(doc['accession'], vdb_strains, adjusted_sequence_accessions)
        if adjusted_sequence_accessions > 0:
            print(str(adjusted_sequence_accessions) + " out of " + str(len(sequences)) + " sequence accessions were adjusted to a match a sequence in " + database + "." + table)

    def adjust_name(self, name, vdb_strains, number_matched):
        '''
        Change strain name if matched to a relaxed strain name in vdb strains
        Return adjusted strain name and total count of matched measurements
        '''
        relax_strain = self.relax_name(name)
        if relax_strain in vdb_strains and vdb_strains[relax_strain] != name:
            number_matched += 1
            return vdb_strains[relax_strain], number_matched
        return name, number_matched

    def link_viruses_to_sequences(self, viruses, sequences):
        '''
        Link the sequence information virus isolate information via the strain name
        '''
        strain_name_to_virus_doc = {}
        for virus in viruses:
            if virus['strain'] not in strain_name_to_virus_doc:
                strain_name_to_virus_doc[virus['strain']] = [virus]
            else:
                strain_name_to_virus_doc[virus['strain']].append(virus)
        for sequence_doc in sequences:
            if sequence_doc['strain'] in strain_name_to_virus_doc:  # determine if sequence has a corresponding virus to link to
                for virus_doc in strain_name_to_virus_doc[sequence_doc['strain']]:
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

    def upload_documents(self, table, documents, database, replace=False, **kwargs):
        if replace:
            print("Deleting documents in database:" + database + "." + table)
            r.table(table).delete().run()
        print("Inserting ", len(documents), "documents")
        self.upload_to_rethinkdb(database, table, documents, **kwargs)

    def upload_to_rethinkdb(self, database, table, documents, overwrite=False, optimal_upload=200, **kwargs):
        if len(documents) > optimal_upload:
            list_documents = [documents[x:x+optimal_upload] for x in range(0, len(documents), optimal_upload)]
        else:
            list_documents = [documents]
        print("Uploading to rethinkdb in " + str(len(list_documents)) + " batches of " + str(optimal_upload) + " documents at a time")
        inserted = 0
        replaced = 0
        for list_docs in list_documents:
            try:
                if not overwrite:
                    document_changes = r.table(table).insert(list_docs, conflict=lambda id, old_doc, new_doc: rethinkdb_updater(id, old_doc, new_doc), return_changes=True).run()
                else:
                    document_changes = r.table(table).insert(list_docs, conflict=lambda id, old_doc, new_doc: rethinkdb_updater_overwrite(id, old_doc, new_doc), return_changes=True).run()
            except:
                raise Exception("Couldn't insert new documents into database", database + "." + table)
            else:
                self.update_timestamp(table, document_changes, **kwargs)
                if document_changes['errors']>0:
                    print("Errors were made when inserting the documents", document_changes['errors'])
                    print(document_changes['errors'])
                inserted += document_changes['inserted']
                replaced += document_changes['replaced']
        print("Ended up inserting " + str(inserted) + " documents into " + database + "." + table)
        print("Ended up updating " + str(replaced) + " documents in " + database + "." + table)

    def update_timestamp(self, table, document_changes, index, **kwargs):
        '''
        Update the timestamp field in the rethink table if changes have been made to the documents
        '''
        updated_documents = []
        if 'changes' in document_changes:
            for doc in document_changes['changes']:
                if doc['new_val'] is not None:
                    updated_documents.append({index: doc['new_val'][index], 'timestamp': self.rethink_io.get_upload_timestamp()})
        if len(updated_documents) > 0:
            r.table(table).insert(updated_documents, conflict='update').run()

    def relaxed_keys(self, indexes, relax_name_method):
        '''
        Use list of indexes
        Create dictionary from relaxed vdb strain names to actual vdb strain names.
        '''
        strains = {}
        for index in indexes:
            strains[relax_name_method(index)] = index
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

def rethinkdb_updater(id, old_doc, new_doc):
    return (new_doc.keys().set_union(old_doc.keys()).map(lambda key:
        r.branch(old_doc.keys().contains(key).and_(new_doc.keys().contains(key).not_()),
            [key, old_doc[key]],
            new_doc.keys().contains(key).and_(old_doc.keys().contains(key).not_()),
            [key, new_doc[key]],
            r.branch(key.eq('sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences'])],
                key.eq('number_sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences']).count()],
                key.eq('timestamp').or_(key.eq('virus_inclusion_date')).or_(key.eq('sequence_inclusion_date')),
                [key, old_doc[key]],
                r.branch(old_doc[key].eq(None).and_(new_doc[key].eq(None).not_()),
                    [key, new_doc[key]],
                    [key, old_doc[key]])
            )
        )
    )).coerce_to('object')

def rethinkdb_updater_overwrite(id, old_doc, new_doc):
    return (new_doc.keys().set_union(old_doc.keys()).map(lambda key:
        r.branch(old_doc.keys().contains(key).and_(new_doc.keys().contains(key).not_()),
            [key, old_doc[key]],
            new_doc.keys().contains(key).and_(old_doc.keys().contains(key).not_()),
            [key, new_doc[key]],
            r.branch(key.eq('sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences'])],
                key.eq('number_sequences'),
                [key, old_doc['sequences'].set_union(new_doc['sequences']).count()],
                key.eq('timestamp').or_(key.eq('virus_inclusion_date')).or_(key.eq('sequence_inclusion_date')),
                [key, old_doc[key]],
                [key, new_doc[key]]
            )
        )
    )).coerce_to('object')
