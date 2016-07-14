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

class gisaid_flu_upload(upload):
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

    def upload(self, preview=False, **kwargs):
        '''
        format virus information, then upload to database
        '''
        print("Uploading Viruses to VDB")
        self.parse(**kwargs)
        self.format()
        self.format_schema()
        self.link_sequences(self.viruses, self.sequences)
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

    def parse(self, path, fname, **kwargs):
        '''
        '''
        fasta_fname = path + fname + ".fasta"
        xls_fname = path + fname + ".xls"
        self.sequences = self.parse_fasta_file(fasta_fname)
        self.viruses = self.parse_gisaid_xls_file(xls_fname, **kwargs)

    def parse_fasta_file(self, fasta, **kwargs):
        '''
        Parse FASTA file with default header formatting
        :return: list of documents(dictionaries of attributes) to upload
        '''
        viruses = []
        try:
            handle = open(fasta, 'r')
        except IOError:
            raise Exception(fasta, "not found")
        else:
            for record in SeqIO.parse(handle, "fasta"):
                content = list(map(lambda x: x.strip(), record.description.replace(">", "").split('|')))
                v = {key: content[ii] if ii < len(content) else "" for ii, key in self.fasta_fields.items()}
                v['sequence'] = str(record.seq)
                v['strain'] = self.fix_name(v['strain'])
                viruses.append(v)
            handle.close()
        return viruses

    def parse_gisaid_xls_file(self, xls, **kwargs):
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
            viruses = [self.fix_gisaid_fields(v, **kwargs) for v in viruses]
            viruses = [self.add_other_attributes(v, **kwargs) for v in viruses]
        return viruses

    def fix_gisaid_fields(self, v, xls_fields_wanted, **kwargs):
        '''
        '''
        # some names don't change because need to be parsed later

        v = {new_field: v[old_field] if old_field in v else None for new_field, old_field in xls_fields_wanted}
        v = {field: self.camelcase_to_snakecase(v[field]) if field in ['gender', 'host'] and v[field] is not None else v[field] for field in v}
        if 'lab' in v:
            v['lab'] = v['lab'].replace(' ', '_').replace('-', '_').lower()
        if 'strain' in v:
            v['strain'] = self.fix_name(v['strain'])
        v = self.fix_age(v)
        v = self.determine_group_fields(v)
        return v

    def fix_age(self, v):
        v['age'] = None
        if 'Host_Age' in v:
            if v['Host_Age'] is not None:
                v['age'] = str(int(v['Host_Age']))
            del v['Host_Age']
        if 'Host_Age_Unit' in v:
            if v['Host_Age_Unit'] is not None and v['age'] is not None:
                v['age'] = v['age'] + v['Host_Age_Unit'].strip().lower()
            del v['Host_Age_Unit']
        return v

    def format(self):
        '''
        format virus information in preparation to upload to database table
        Flu needs to also format country
        '''
        print('Formatting for upload')
        self.define_countries()
        self.define_regions()
        for virus in self.viruses:
            self.format_date(virus)
            self.format_country(virus)
            self.format_region(virus)
            self.format_place(virus)

    def fix_name(self, name):
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '')
        split_name = tmp_name.split('/')
        if split_name[1].isupper():
            split_name[1] = split_name[1].title()  # B/WAKAYAMA-C/2/2016 becomes B/Wakayama-C/2/2016
        split_name[2] = split_name[2].lstrip('0')  # A/Mali/013MOP/2015 becomes A/Mali/13MOP/2015
        if '-' in split_name[1]:
            print(name)
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
            if v['host'] == "human":
                v['location'] = v['strain'].split('/')[1]
            elif len(v['strain'].split('/')) ==4:
                v['location'] = v['strain'].split('/')[1]
            elif len(v['strain'].split('/')) == 5:
                v['location'] = v['strain'].split('/')[2]
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
                v['country'] = None
                print("couldn't parse country for", v['strain'])

    def determine_group_fields(self, v, **kwargs):
        '''
        Determine and assign genetic group fields
        '''
        # determine virus type from strain name
        temp_subtype = v['Subtype']
        temp_lineage = v['Lineage']
        del v['Subtype']
        del v['Lineage']
        v['vtype'] = 'undetermined'
        v['subtype'] = 'undetermined'
        v['lineage'] = 'undetermined'
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
                else:
                    print("Couldn't parse virus type from strain name: " + v['strain'], (v['subtype'],v['lineage']))
        return v

    def format_schema(self):
        '''
        move sequence information into nested 'sequences' field
        '''
        for virus in self.viruses:
            virus['sequences'] = []

    def link_sequences(self, viruses, sequences):
        '''
        Link the sequence information from the fasta file to the virus isolate information from the xls file
        '''
        strain_name_to_virus = {virus['strain']: virus for virus in viruses}
        for sequence_info in sequences:
            strain = sequence_info['strain']
            virus = strain_name_to_virus[strain]
            del sequence_info['strain']
            if 'sequences' not in virus:
                virus['sequences'] = []
            virus['sequences'].append(sequence_info)

    def filter(self):
        '''
        Filter out viruses without correct format, check  optional and upload attributes
        '''
        print(str(len(self.viruses)) + " viruses before filtering")
        self.rethink_io.check_optional_attributes(self.viruses, [])
        self.viruses = filter(lambda v: self.rethink_io.check_required_attributes(v, self.upload_fields, self.index_field), self.viruses)
        print(str(len(self.viruses)) + " viruses after filtering")
        self.viruses = filter(lambda v: 'sequences' in v and len(v['sequences']) > 0, self.viruses)
        print(str(len(self.viruses)) + " viruses after filtering")

    def update_document_sequence(self, document, v, **kwargs):
        '''
        Update sequence fields if matching accession or sequence
        Append sequence to sequence list if no matching accession or sequence
        If sequence is currently empty, ie accession: null, sequence: null, then replace
        '''
        updated = False
        if 'sequences' in document and 'sequences' in v:
            doc_seqs = document['sequences']
            virus_seqs = v['sequences']
            for virus_seq in virus_seqs:
                if len(doc_seqs) == 0:
                    doc_seqs.extend(virus_seqs)
                    updated = True
                if virus_seq['accession'] is not None: # try comparing accession's first
                    if all(virus_seq['accession'] != seq_info['accession'] for seq_info in doc_seqs):
                        updated = self.append_new_sequence(document, virus_seq)
                    else:
                        updated = self.update_sequence_citation_field(document, v, 'accession', self.sequence_upload_fields+self.sequence_optional_fields, self.citation_optional_fields, **kwargs)
        return updated

    def append_new_sequence(self,document, virus_seq):
        '''
        New sequence information, so append to document sequences and citations field
        '''
        document['sequences'].append(virus_seq)
        return True

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0: 'strain', 1: 'isolate_id', 2: 'accession', 3:'locus', 4: 'passage'}
    #              >>B/Austria/896531/2016  | EPI_ISL_206054 | 687738 | HA | Siat 1
    setattr(args, 'fasta_fields', fasta_fields)
    xls_fields_wanted = [('strain', 'Isolate_Name'), ('isolate_id', 'Isolate_Id'), ('date', 'Collection_Date'),
                             ('host', 'Host'), ('Subtype', 'Subtype'), ('Lineage', 'Lineage'),
                             ('Location', 'Location'), ('lab', 'Submitting_Lab'), ('Host_Age', 'Host_Age'),
                             ('Host_Age_Unit', 'Host_Age_Unit'), ('gender', 'Host_Gender')]
    setattr(args, 'xls_fields_wanted', xls_fields_wanted)
    if args.path is None:
        args.path = "vdb/data/" + args.virus + "/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = gisaid_flu_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)