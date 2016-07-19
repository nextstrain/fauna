import os, re, time, datetime, csv, sys
import numpy as np
import rethinkdb as r
from Bio import SeqIO
from Bio import AlignIO
from upload import upload
from upload import parser

parser.add_argument('--vtype', default=None, help="type of virus, if applicable")
parser.add_argument('--subtype', default=None, help="subtype of virus")
parser.add_argument('--lineage', default=None, help="lineage of virus")

class flu_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.grouping_upload_fields = ['vtype', 'subtype', 'lineage']

        # patterns from the subtype and lineage fields in the GISAID fasta file
        self.patterns = {('a / h3n2', ''): ('a', 'h3n2', 'seasonal_h3n2'),
                    ('a / h3n0', ''): ('a', 'h3n2', 'seasonal_h3n2'),  # often mistake in GISAID file, aligned closest to H3N2
                    ('a / h1n1', 'pdm09'): ('a', 'h1n1', 'seasonal_h1n1pdm'),
                    ('b / h0n0', 'Victoria'): ('b', 'undetermined', 'seasonal_vic'),
                    ('b / h0n0', 'Yamagata'): ('b', 'undetermined', 'seasonal_yam'),
                    ('a / h1n1', 'seasonal'): ('a', 'h1n1', 'seasonal_h1n1'),
                    ('a / h7n9', ''): ('a', 'h7n9', 'undetermined'),
                    ('a / h5n1', ''): ('a', 'h5n1', 'undetermined'),
                    ('a / h6n1', ''): ('a', 'h6n1', 'undetermined'),
                    ('a / h5n6', ''): ('a', 'h5n6', 'undetermined')}
        self.outgroups = {lineage: SeqIO.read('source-data/'+lineage+'_outgroup.gb', 'genbank') for lineage in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']}
        self.outgroup_patterns = {'H3N2': ('a', 'h3n2', 'seasonal_h3n2'),
                                  'H1N1pdm': ('a', 'h1n1', 'seasonal_h1n1pdm'),
                                  'Vic': ('b', 'undetermined', 'seasonal_vic'),
                                  'Yam': ('b', 'undetermined', 'seasonal_yam')}
        db_viruses = r.db(self.database).table(self.virus).run()
        self.db_strains = {v['strain'] for v in db_viruses}
        self.virus_to_sequence_transfer_fields = ['submission_date']

    def parse(self, path, fname, **kwargs):
        '''
        '''
        fasta_fname = path + fname + ".fasta"
        xls_fname = path + fname + ".xls"
        sequences = self.parse_fasta_file(fasta_fname, **kwargs)
        viruses = self.parse_gisaid_xls_file(xls_fname, **kwargs)
        print("Parsed " + str(len(viruses)) + " viruses and " + str(len(sequences)) + " sequences from file")
        return (viruses, sequences)

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
                v = {key: content[ii] if ii < len(content) else "" for ii, key in sequence_fasta_fields.items()}
                v['sequence'] = str(record.seq)
                v['strain'] = self.fix_name(v['strain'])
                v = self.add_sequence_fields(v, **kwargs)
                sequences.append(v)
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
        return viruses

    def format(self, documents, **kwargs):
        '''
        format virus information in preparation to upload to database table
        Flu needs to also format country
        '''
        self.define_countries()
        self.define_regions()
        for doc in documents:
            if 'strain' in doc:
                doc['strain'] = self.fix_name(doc['strain'])
            self.fix_casing(doc)
            self.fix_age(doc)
            self.add_virus_fields(doc, **kwargs)
            self.determine_group_fields(doc, **kwargs)
            self.format_date(doc)
            self.format_country(doc)
            self.format_region(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.format_place(doc)

    def fix_casing(self, doc):
        for field in ['originating_lab', 'submitting_lab']:
            if field in doc and doc[field] is not None:
                doc[field] = doc[field].replace(' ', '_').replace('-', '_').lower()
        for field in ['gender', 'host', 'locus']:
            if field in doc and doc[field] is not None:
                doc[field] = self.camelcase_to_snakecase(doc[field])
        if 'accession' in doc and doc['accession'] is not None:
            doc['accession'] = 'EPI' + doc['accession']

    def fix_age(self, doc):
        doc['age'] = None
        if 'Host_Age' in doc:
            if doc['Host_Age'] is not None:
                doc['age'] = str(int(doc['Host_Age']))
            del doc['Host_Age']
        if 'Host_Age_Unit' in doc:
            if doc['Host_Age_Unit'] is not None and doc['age'] is not None:
                doc['age'] = doc['age'] + doc['Host_Age_Unit'].strip().lower()
            del doc['Host_Age_Unit']
        return doc

    def fix_name(self, name):
        if '(' in name and ')' in name:  # A/Eskisehir/359/2016 (109) -> A/Eskisehir/359/2016 ; A/South Australia/55/2014  IVR145  (14/232) -> A/South Australia/55/2014  IVR145
            name = re.match(r'^([^(]+)', name).group(1)
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '').replace('&', '')
        split_name = tmp_name.split('/')
        if split_name[1].isupper():
            split_name[1] = split_name[1].title()  # B/WAKAYAMA-C/2/2016 becomes B/Wakayama-C/2/2016
        split_name[2] = split_name[2].lstrip('0')  # A/Mali/013MOP/2015 becomes A/Mali/13MOP/2015
        result_name = '/'.join(split_name)
        return result_name

    def define_countries(self):
        '''
        open synonym to country dictionary
        Location is to the level of country of administrative division when available
        '''
        reader = csv.DictReader(open("source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
        self.label_to_country = {}
        for line in reader:
            self.label_to_country[line['label'].lower()] = line['country']

    def format_country(self, v):
        '''
        Label viruses with country based on strain name
        '''
        if "country" not in v or v['country'].lower() == 'unknown':
            v['country'] = '?'
        if len(v['strain'].split('/')) ==4:
            v['location'] = v['strain'].split('/')[1]
        elif len(v['strain'].split('/')) == 5:
            v['location'] = v['strain'].split('/')[2]
        else:
            print("couldn't parse country for", v['strain'], v['gisaid_location'])
            v['location'] = None
        try:
            label = re.match(r'^([^/]+)', v['location']).group(1).lower()						# check first for whole geo match
            if label in self.label_to_country:
                v['country'] = self.label_to_country[label]
            else:
                label = re.match(r'^([^\-^\/]+)', v['location']).group(1).lower()			# check for partial geo match A/CHIBA-C/61/2014
            if label in self.label_to_country:
                v['country'] = self.label_to_country[label]
            else:
                label = re.match(r'^([A-Z][a-z]+)[A-Z0-9]', v['location']).group(1).lower()			# check for partial geo match
            if label in self.label_to_country:
                v['country'] = self.label_to_country[label]
        except:
            if 'gisaid_location' in v:
                if len(v['gisaid_location'].split("/")) > 1:
                    v['country'] = v['gisaid_location'].split("/")[1].strip()
                else:
                    print("couldn't parse country for", v['strain'], v['gisaid_location'])

    def determine_group_fields(self, v, **kwargs):
        '''
        Determine and assign genetic group fields
        '''
        # determine virus type from strain name
        if 'Subtype' in v and 'Lineage' in v:
            if v['Subtype'] is not None:
                temp_subtype = v['Subtype']
            else:
                temp_subtype = ''
            del v['Subtype']
            if v['Lineage'] is not None:
                temp_lineage = v['Lineage']
            else:
                temp_lineage = ''
            del v['Lineage']
            v['vtype'] = 'tbd'
            v['subtype'] = 'tbd'
            v['lineage'] = 'tbd'
            if (temp_subtype, temp_lineage) in self.patterns:  #look for pattern from GISAID fasta file
                match = self.patterns[(temp_subtype, temp_lineage)]
                v['vtype'] = match[0].lower()
                v['subtype'] = match[1].lower()
                v['lineage'] = match[2].lower()
            else:
                if '/' in v['strain']:
                    #determine virus type from strain name
                    fields = v['strain'].split('/')
                    if re.match(r'^[abcd]$',fields[0].strip()):
                        v['vtype'] = fields[0].strip().lower()
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
        args.path = "vdb/data/" + args.virus + "/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = flu_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)