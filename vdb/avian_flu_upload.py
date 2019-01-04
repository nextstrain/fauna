import os, re, time, datetime, csv, sys, json
import numpy as np
import rethinkdb as r
from Bio import SeqIO
from Bio import AlignIO
from upload import upload
from upload import get_parser
from unidecode import unidecode

parser = get_parser()
parser.add_argument('--upload_directory', default=False, action="store_true", help='upload all xls and fasta files in directory')
parser.add_argument('--vtype', default=None, help="type of virus, if applicable")
parser.add_argument('--subtype', default=None, help="subtype of virus")
parser.add_argument('--lineage', default=None, help="lineage of virus")

class flu_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.grouping_upload_fields = ['vtype', 'subtype', 'lineage']
        # patterns from the subtype and lineage fields in the GISAID fasta file
        self.patterns = {('a / h1n1', 'pdm09'): ('a', 'h1n1', 'seasonal_h1n1pdm'),
                    ('a / h1n2', ''): ('a', 'h1n2', None),
                    ('a / h1n2', 'seasonal'): ('a', 'h1n2', 'seasonal_h1n2'),
                    ('a / h2n2', ''): ('a', 'h2n2', None),
                    ('a / h3n2', ''): ('a', 'h3n2', 'seasonal_h3n2'),
                    ('a / h3n2', 'seasonal'): ('a', 'h3n2', 'seasonal_h3n2'),
                    ('a / h3n3', ''): ('a', 'h3n3', None),
                    ('a / h5n1', ''): ('a', 'h5n1', None),
                    ('a / h5n6', ''): ('a', 'h5n6', None),
                    ('a / h6n1', ''): ('a', 'h6n1', None),
                    ('a / h7n1', ''): ('a', 'h7n1', None),
                    ('a / h7n2', ''): ('a', 'h7n2', None),
                    ('a / h7n3', ''): ('a', 'h7n3', None),
                    ('a / h7n7', ''): ('a', 'h7n7', None),
                    ('a / h7n9', ''): ('a', 'h7n9', None),
                    ('a / h9n2', ''): ('a', 'h9n2', None),
                    ('a / h10n7', ''): ('a', 'h10n7', None),
                    ('a / h10n8', ''): ('a', 'h10n8', None),
                    ('a / h11', ''): ('a', 'h11', None),
                    ('b / h0n0', 'victoria'): ('b', None, 'seasonal_vic'),
                    ('b / h0n0', 'yamagata'): ('b', None, 'seasonal_yam'),
                    ('b', 'victoria'): ('b', None, 'seasonal_vic'),
                    ('b', 'yamagata'): ('b', None, 'seasonal_yam')}
        self.outgroups = {lineage: SeqIO.read('source-data/'+lineage+'_outgroup.gb', 'genbank') for lineage in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']}
        self.outgroup_patterns = {'H3N2': ('a', 'h3n2', 'seasonal_h3n2'),
                                  'H1N1': ('a', 'h1n1', 'seasonal_h1n1'),
                                  'H1N1pdm': ('a', 'h1n1', 'seasonal_h1n1pdm'),
                                  'Vic': ('b', None, 'seasonal_vic'),
                                  'Yam': ('b', None, 'seasonal_yam')}
        self.strain_fix_fname = "source-data/flu_strain_name_fix.tsv"
        self.location_fix_fname = "source-data/flu_location_fix.tsv"
        self.virus_to_sequence_transfer_fields = ['submission_date']
        self.fix = set()

    def parse(self, path, fname, upload_directory, **kwargs):
        '''
        '''
        viruses, sequences = [], []
        if upload_directory:
            import glob
            for xls_fname, fasta_fname in zip(glob.glob(path + "gisaid*.xls"), glob.glob(path + "gisaid*.fasta")):
                parsed = self.parse_files(xls_fname, fasta_fname, **kwargs)
                viruses.extend(parsed[0])
                sequences.extend(parsed[1])
        else:
            fasta_fname = path + fname + ".fasta"
            xls_fname = path + fname + ".xls"
            viruses, sequences = self.parse_files(xls_fname, fasta_fname, **kwargs)
        print("Parsed total of " + str(len(viruses)) + " viruses and " + str(len(sequences)) + " sequences from files")
        return viruses, sequences

    def parse_files(self, xls_fname, fasta_fname, **kwargs):
        '''
        parse linked xls and fasta downloaded from gisaid
        '''
        viruses = self.parse_gisaid_xls_file(xls_fname, **kwargs)
        sequences = self.parse_fasta_file(fasta_fname, **kwargs)
        print("Parsed " + str(len(viruses)) + " viruses and " + str(len(sequences)) + " sequences from files", fasta_fname, xls_fname)
        return viruses, sequences

    def parse_fasta_file(self, fasta, **kwargs):
        '''
        Parse FASTA file with default header formatting
        :return: list of documents(dictionaries of attributes) to upload
        '''
        sequences = []
        try:
            handle = open(fasta, 'r')
        except IOError:
            raise Exception(fasta, "not found")
        else:
            for record in SeqIO.parse(handle, "fasta"):
                content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
                s = {key: content[ii] if ii < len(content) else "" for ii, key in sequence_fasta_fields.items()}
                s['sequence'] = str(record.seq)
                s = self.add_sequence_fields(s, **kwargs)
                sequences.append(s)
            handle.close()
        return sequences

    def parse_gisaid_xls_file(self, xls, xls_fields_wanted, **kwargs):
        '''
        parse excel file using pandas
        :return: list of documents(dictionaries of attributes) to upload
        '''
        import pandas
        try:
            handle = open(xls, 'rb')
        except IOError:
            raise Exception(xls, "not found")
        else:
            df = pandas.read_excel(handle)
            df = df.where((pandas.notnull(df)), None)  # convert Nan type to None
            viruses = df.to_dict('records')
            viruses = [{new_field: v[old_field] if old_field in v else None for new_field, old_field in xls_fields_wanted} for v in viruses]
            viruses = [self.add_virus_fields(v, **kwargs) for v in viruses]
        return viruses

    def format_viruses(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        if self.strain_fix_fname is not None:
            self.fix_whole_name = self.define_strain_fixes(self.strain_fix_fname)
        if self.location_fix_fname is not None:
            self.fix_location = self.define_location_fixes(self.location_fix_fname)
        self.define_countries("source-data/geo_synonyms.tsv")
        self.define_regions("source-data/geo_regions.tsv")
        self.define_location_label_fixes("source-data/flu_fix_location_label.tsv")
        for doc in documents:
            if 'strain' in doc:
                doc['strain'], doc['gisaid_strain'] = self.fix_name(doc['strain'])
            else:
                print("Missing strain name!")
            self.fix_casing(doc)
            self.fix_age(doc)
            self.format_host(doc)
            self.determine_group_fields(doc, self.patterns)
            self.format_date(doc)
            self.format_country(doc) # first format from strain name
            if self.fix_location is not None: # override with fixes
                if doc['strain'] in self.fix_location:
                    doc['location'] = self.fix_location[doc['strain']]
            self.format_place(doc, determine_location=True)
            self.format_region(doc)
            self.rethink_io.check_optional_attributes(doc, [])

    def format_sequences(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        for doc in documents:
            if 'strain' in doc:
                doc['strain'], doc['gisaid_strain'] = self.fix_name(doc['strain'])
            else:
                print("Missing strain name!")
            self.format_date(doc)
            self.format_passage(doc, 'passage', 'passage_category')
            self.format_passage(doc, 'virus_strain_passage', 'virus_strain_passage_category') #BP
            self.format_passage(doc, 'serum_antigen_passage', 'serum_antigen_passage_category') #BP
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)
        print("Names that need to be fixed")
        for name in sorted(self.fix):
            print(name)

    def filter(self, documents, index, **kwargs):
        '''
        filter out certain documents
        '''
        print(str(len(documents)) + " documents before filtering")
        documents = filter(lambda doc: index in doc, documents)
        remove_labels = []
        # remove certain documents from gisaid files that were not actually isolated from humans
        result_documents = [doc for doc in documents if all(label not in doc['strain'] for label in remove_labels)]
        #result_documents = [doc for doc in result_documents if self.correct_strain_format(doc['strain'], doc['gisaid_strain'])]
        print(str(len(result_documents)) + " documents after filtering")
        return result_documents

    def correct_strain_format(self, strain, original_strain):
        # Okay Patterns: B/Brisbane/46/2015, A/HongKong/1968, A/Zambia/13/176/2013 or A/Cologne/Germany/12/2009 or A/Algeria/G0164/15/2015 or A/India/Delhi/DB106/2009, A/Cameroon/LEID/01/11/1387/2011, A/India/M/Enc/1/2003
        if re.match(r'[A|B]/[A-Za-z-]+/([A-Za-z0-9_-]+/)*[0-9]{4}$', strain) or re.match(r'[A|B]/[A-Za-z-]+/([A-Za-z0-9_-]+/){2}[0-9]{4}$', strain)\
                or re.match(r'[A|B]/([A-Za-z-]+/){2}([0-9]+/){3}[0-9]{4}$', strain):
            return True
        else:
            print("This strain name was not in the correct format and will be filtered out", strain, original_strain)
            self.fix.add(strain)

    def fix_casing(self, doc):
        '''
        fix gisaid specific fields casing
        '''
        for field in ['originating_lab', 'submitting_lab']:
            if field in doc and doc[field] is not None:
                doc[field] = doc[field].replace(' ', '_').replace('-', '_').lower()
        for field in ['gender', 'host', 'locus']:
            if field in doc and doc[field] is not None:
                doc[field] = self.camelcase_to_snakecase(doc[field])
        if 'accession' in doc and doc['accession'] is not None:
            doc['accession'] = 'EPI' + doc['accession']

    def fix_age(self, doc):
        '''
        Combine gisaid age information into one age field
        '''
        temp_age, temp_age_unit = None, None
        doc['age'] = None
        if 'Host_Age' in doc:
            try:
                temp_age = str(int(float(doc['Host_Age'])))
            except:
                pass
            del doc['Host_Age']
        if 'Host_Age_Unit' in doc:
            if isinstance(doc['Host_Age_Unit'], basestring):
                temp_age_unit = doc['Host_Age_Unit'].lower()
            else:
                temp_age_unit = 'y'
            del doc['Host_Age_Unit']
        if isinstance(temp_age, basestring) and isinstance(temp_age_unit, basestring):
            doc['age'] = temp_age + temp_age_unit
        return doc

    def define_location_fixes(self, fname):
        '''
        Open location fix file and define corresponding dictionaries
        '''
        reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
        fix_location = {}
        for line in reader:
            fix_location[line['label'].decode('unicode-escape')] = line['fix']
        return fix_location

    def define_location_label_fixes(self, fname):
        reader = csv.DictReader(filter(lambda row: row[0]!='#', open(fname)), delimiter='\t')
        self.label_to_fix = {}
        for line in reader:
            self.label_to_fix[line['label'].decode('unicode-escape').replace(' ', '').lower()] = line['fix']

    def fix_name(self, name):
        '''
        Fix strain names
        '''
        # replace all accents with ? mark
        original_name = name.encode('ascii', 'replace')
        # return original_name, original_name
        # Replace whole strain names
        name = self.replace_strain_name(original_name, self.fix_whole_name)
        name = name.replace('H1N1', '').replace('H5N6', '').replace('H3N2', '').replace('H5N1', '').replace('H7N9', '')\
            .replace('Influenza A Virus', '').replace('segment 4 hemagglutinin (HA) gene', '').replace("segment 6 neuraminidase (NA) gene", "")\
            .replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '').replace('&', '').replace(' ', '')\
            .replace('\'', '').replace('>', '').replace('-like', '').replace('+', '')
        name = name.lstrip('-').lstrip('_').lstrip(')').lstrip('(')
        name = name.lstrip('-').rstrip('_').rstrip(')').rstrip('(')

        split_name = name.split('/')
        # check location labels in strain names for fixing
        for index, label in enumerate(split_name):
            if label.replace(' ', '').lower() in self.label_to_fix:
                split_name[index] = self.label_to_fix[label.replace(' ', '').lower()]
        name = '/'.join(split_name)
        name = self.flu_fix_patterns(name)

        # Strip leading zeroes, change all capitalization location field to title case
        split_name = name.split('/')
        if len(split_name) == 4:
            if split_name[1].isupper() or split_name[1].islower():
                split_name[1] = split_name[1].title()  # B/WAKAYAMA-C/2/2016 becomes B/Wakayama-C/2/2016
            split_name[2] = split_name[2].lstrip('0')  # A/Mali/013MOP/2015 becomes A/Mali/13MOP/2015
            split_name[3] = split_name[3].lstrip('0')  # A/Cologne/Germany/01/2009 becomes A/Cologne/Germany/1/2009
        result_name = '/'.join(split_name).strip()

        return result_name, original_name

    def flu_fix_patterns(self, name):
        # various name patterns that need to be fixed
        # capitalization of virus type
        if re.match(r'([a|b])([\w\s\-/]+)', name):  #b/sydney/508/2008    B/sydney/508/2008
            name = re.match(r'([a|b])([\w\s\-/]+)', name).group(1).upper() + re.match(r'([a|b])([\w\s\-/]+)', name).group(2)
        # remove inner parentheses and their contents
        if re.match(r'([^(]+)[^)]+\)(.+)', name):  # A/Egypt/51(S)/2006
            name = re.match(r'([^(]+)[^)]+\)(.+)', name).group(1) + re.match(r'([^(]+)[^)]+\)(.+)', name).group(2)
        # remove ending parentheses and their contents
        if re.match(r'([^(]+)[^)]+\)$', name):  # A/Eskisehir/359/2016 (109) -> A/Eskisehir/359/2016 ; A/South Australia/55/2014  IVR145  (14/232) -> A/South Australia/55/2014  IVR145
            name = re.match(r'([^(]+)[^)]+\)$', name).group(1)
        # Strip trailing slashes
        name = name.rstrip('/')  # A/NorthernTerritory/60/68//  A/Paris/455/2015/
        # Change two digit years to four digit years
        if re.match(r'([\w\s\-/]+)/([0-9][0-9])$', name):  #B/Florida/1/96 -> B/Florida/1/1996
            year = re.match(r'([\w\s\-/]+)/([0-9][0-9])$', name).group(2)
            if int(year) < 66:
                name = re.match(r'([\w\s\-/]+)/([0-9][0-9])$', name).group(1) + "/20" + year
            else:
                name = re.match(r'([\w\s\-/]+)/([0-9][0-9])$', name).group(1) + "/19" + year
        return name

    def format_host(self, v):
        '''
        Fix host formatting
        '''
        if v['host'] is not None:
            if v['host'] == "chicken":
                v['host'] = "avian"
            if v['host'] == "duck":
                v['host'] = "avian"
            if v['host'] == "gallusgallus":
                v['host'] = "avian"
            if v['host'] == "gallusgallusdomesticus":
                v['host'] = "avian"
            if v['host'] == "anasclypeata":
                v['host'] = "avian"
            if v['host'] == "anasplatyrhynchos":
                v['host'] = "avian"
            if v['host'] == "anassp.":
                v['host'] = "avian"
            if v['host'] == "anaspoecilorhyncha":
                v['host'] = "avian"
            if v['host'] == "anasdiscors":
                v['host'] = "avian"
            if v['host'] == "anascrecca":
                v['host'] = "avian"
            if v['host'] == "anasstrepera":
                v['host'] = "avian"
            if v['host'] == "passerine":
                v['host'] = "avian"
            if v['host'] == "larusridibundus":
                v['host'] = "avian"
            if v['host'] == "anascarolinensis":
                v['host'] = "avian"
            if v['host'] == "us_quail":
                v['host'] = "avian"
            if v['host'] == "goose":
                v['host'] = "avian"
            if v['host'] == "anasrubripes":
                v['host'] = "avian"
            if v['host'] == "anasamericana":
                v['host'] = "avian"
            if v['host'] == "corvus":
                v['host'] = "avian"
            if v['host'] == "falcoperegrinus":
                v['host'] = "avian"
            if v['host'] == "zosteropsjaponicus":
                v['host'] = "avian"
            if v['host'] == "cygnuscygnus":
                v['host'] = "avian"
            if v['host'] == "falcon":
                v['host'] = "avian"
            if v['host'] == "eagle":
                v['host'] = "avian"
            if v['host'] == "turkey":
                v['host'] = "avian"
            if v['host'] == "graculareligiosa":
                v['host'] = "avian"
            if v['host'] == "chencanagica":
                v['host'] = "avian"
            if v['host'] == "anserindicus":
                v['host'] = "avian"
            if v['host'] == "passermontanus":
                v['host'] = "avian"
            if v['host'] == "arenariainterpres":
                v['host'] = "avian"
            if v['host'] == "otheravian":
                v['host'] = "avian"
            if v['host'] == "avian":
                v['host'] = "avian"
            if v['host'] == "coturnix":
                v['host'] = "avian"
            if v['host'] == "guineafowl":
                v['host'] = "avian"
            if v['host'] == "cairinamoschata":
                v['host'] = "avian"
            if v['host'] == "anascyanoptera":
                v['host'] = "avian"
            if v['host'] == "feline":
                v['host'] = "nonhuman_mammal"
            if v['host'] == "watersample":
                v['host'] = "environment"
            if v['host'] == "surfaceswab":
                v['host'] = "environment"
            if v['host'] == "feces":
                v['host'] = "environment"
            if v['host'] == "watersample":
                v['host'] = "environment"

    def format_country(self, v):
        '''
        Label viruses with country based on strain name
        A/Taiwan/1/2013 is human virus. Four fields total. Take second field.
        A/Chicken/Taiwan/1/2013 is animal virus. Five field total. Take third field.
        Else, take GISAID location.
        '''
        strain_name = v['strain']
        original_name = v['gisaid_strain']
        result = None
        field_count = 0
        if '/' in strain_name:
            field_count = len(strain_name.split('/'))
        if field_count == 4:
            loc = strain_name.split('/')[1].replace(" ", "")
            result = self.determine_location(loc)
        elif field_count == 5:
            loc = strain_name.split('/')[2].replace(" ", "")
            result = self.determine_location(loc)

        if v['gisaid_location'] is not None and result is None:
            loc = v['gisaid_location'].split('/')[-1].replace(" ", "")
            result = self.determine_location(loc)

        if result is not None:
            v['location'], v['division'], v['country'] = result
        else:
            v['location'], v['division'], v['country'] = None, None, None
            print("couldn't parse country for ", strain_name, "gisaid location", v['gisaid_location'], original_name)

        if v['division'] == v['country']:
            v['division'] == '?'

    def format_passage(self, doc, initial_field, new_field, **kwargs):
        '''
        Separate passage into general categories
        Regex borrowed from McWhite et al. 2016
        '''
        if initial_field in doc and doc[initial_field] is not None:
            passage = doc[initial_field].upper()
            passage_category = "undetermined"
            if re.search(r'AM[1-9]|E[1-9]|AMNIOTIC|EGG|EX|AM_[1-9]', passage):   # McWhite
                passage_category = "egg"
            elif re.search(r'AM-[1-9]|EMBRYO|^E$', passage):
                passage_category = "egg"
            elif re.search(r'LUNG|P0|OR_|ORIGINAL|CLINICAL|DIRECT', passage):    # McWhite
                passage_category = "unpassaged"
            elif re.search(r'ORGINAL|ORIGNAL|CLINCAL|THROAT|PRIMARY|NASO|AUTOPSY|BRONCHIAL|INITIAL|NASAL|NOSE|ORIG|SWAB', passage):
                passage_category = "unpassaged"
            elif re.search(r'TMK|RMK|RHMK|RII|PMK|R[1-9]|RX', passage):    # McWhite
                passage_category = "cell"
            elif re.search(r'S[1-9]|SX|SIAT|MDCK|MCDK|C[1-9]|CX|M[1-9]|MX|X[1-9]|^X_$', passage):  # McWhite
                passage_category = "cell"
            elif re.search(r'C_[1-9]|C [1-9]|MD[1-9]|MK[1-9]|MEK[1-9]', passage):
                passage_category = "cell"
            elif re.search(r'[Cc][Ee][Ll][Ll]', passage):
                passage_category = "cell"
            elif re.search(r'^S[1-9]_$| ^SX_$|SIAT2_SIAT1|SIAT3_SIAT1', passage):    # McWhite
                passage_category = "cell"
            elif re.search(r'UNKNOWN|UNDEFINED|NOT SPECIFIED|DIFFERENT ISOLATION SOURCES', passage):
                pass
            doc[new_field] = passage_category
        else:
            doc[initial_field] = None
            doc[new_field] = None

    def determine_group_fields(self, v, patterns, **kwargs):
        '''
        Determine and assign genetic group fields
        '''
        # determine virus type from strain name
        v['vtype'], v['subtype'], v['lineage'] = 'tbd', 'tbd', 'tbd'
        temp_subtype = ''
        temp_lineage = ''
        if 'Subtype' in v:
            if v['Subtype'] is not None:
                temp_subtype = v['Subtype'].lower()
            del v['Subtype']
        if 'Lineage' in v:
            if v['Lineage'] is not None:
                temp_lineage = v['Lineage'].lower()
            del v['Lineage']
        if (temp_subtype, temp_lineage) in patterns:  #look for pattern from GISAID fasta file
                match = patterns[(temp_subtype, temp_lineage)]
                v['vtype'], v['subtype'], v['lineage'] = match[0], match[1], match[2]
        return v

if __name__=="__main__":
    args = parser.parse_args()
    sequence_fasta_fields = {0: 'accession', 1: 'strain', 2: 'isolate_id', 3:'locus', 4: 'passage', 5: 'submitting_lab'}
    #              >>B/Austria/896531/2016  | EPI_ISL_206054 | 687738 | HA | Siat 1
    setattr(args, 'fasta_fields', sequence_fasta_fields)
    xls_fields_wanted = [('strain', 'Isolate_Name'), ('isolate_id', 'Isolate_Id'), ('collection_date', 'Collection_Date'),
                             ('host', 'Host'), ('Subtype', 'Subtype'), ('Lineage', 'Lineage'),
                             ('gisaid_location', 'Location'), ('originating_lab', 'Originating_Lab'), ('Host_Age', 'Host_Age'),
                             ('Host_Age_Unit', 'Host_Age_Unit'), ('gender', 'Host_Gender'), ('submission_date', 'Submission_Date')]
    setattr(args, 'xls_fields_wanted', xls_fields_wanted)
    if args.path is None:
        args.path = "data/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = flu_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
