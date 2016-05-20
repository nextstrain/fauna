import os, re, time, datetime, csv, sys
import numpy as np
import rethinkdb as r
from Bio import SeqIO
from Bio import AlignIO
from vdb_upload import vdb_upload
from vdb_upload import parser

parser.add_argument('--vtype', default=None, help="type of virus, if applicable")
parser.add_argument('--subtype', default=None, help="subtype of virus")
parser.add_argument('--lineage', default=None, help="lineage of virus")
class Flu_vdb_upload(vdb_upload):

    def __init__(self, **kwargs):
        vdb_upload.__init__(self, **kwargs)
        self.grouping_upload_fields = ['vtype', 'subtype', 'lineage']

        # patterns from the subtype and lineage fields in the GISAID fasta file
        self.patterns = {('a / h3n2', ''): ('a', 'h3n2', 'seasonal_h3n2'),
                    ('a / h3n0', ''): ('a', 'h3n2', ''), # often mistake in GISAID file, aligned closest to H3N2
                    ('a / h1n1', 'pdm09'): ('a', 'h1n1', 'seasonal_h1n1pdm'),
                    ('b / h0n0', 'Victoria'): ('b', 'undetermined', 'seasonal_vic'),
                    ('b / h0n0', 'Yamagata'): ('b', 'undetermined', 'seasonal_yam'),
                    ('a / h1n1', 'seasonal'): ('a', 'h1n1', 'seasonal_h1n1'),
                    ('a / h7n9', ''): ('a', 'h7n9', 'undetermined'),
                    ('a / h5n1', ''): ('a', 'h5n1', 'undetermined'),
                    ('a / h6n1', ''): ('a', 'h6n1', 'undetermined'),
                    ('a / h5n6', ''): ('a', 'h5n6', 'undetermined')}
        self.outgroups = {lineage: SeqIO.read('vdb/source-data/'+lineage+'_outgroup.gb', 'genbank') for lineage in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']}
        self.outgroup_patterns = {'H3N2': ('a', 'h3n2', ''),
                                  'H1N1pdm': ('a', 'h1n1', 'pdm09'),
                                  'Vic': ('b', 'undetermined', 'victoria'),
                                  'Yam': ('b', 'undetermined', 'yamagata')}
        db_viruses = r.db(self.database).table(self.virus).run()
        self.db_strains = {v['strain'] for v in db_viruses}

    def parse_fasta_file(self, fasta, **kwargs):
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
                v['sequence'] = str(record.seq)
                self.add_other_attributes(v, **kwargs)
                self.determine_group_fields(v, record, **kwargs)
                viruses.append(v)
            handle.close()
        return viruses

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
            self.format_place(virus)
            self.format_region(virus)
            self.delete_extra_fields(virus)

    def define_countries(self):
        '''
        open synonym to country dictionary
        Location is to the level of country of administrative division when available
        '''
        reader = csv.DictReader(open("vdb/source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
        self.label_to_country = {}
        for line in reader:
            self.label_to_country[line['label'].lower()] = line['country'].lower()

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

    def determine_group_fields(self, v, seq, vtype, subtype, lineage, **kwargs):
        '''
        Determine and assign genetic group fields
        '''
        # determine virus type from strain name
        temp_subtype = v['subtype'].replace('_', ' ')
        temp_lineage = v['lineage'].replace('_', ' ')
        v['vtype'] = 'undetermined'
        v['subtype'] = 'undetermined'
        v['lineage'] = 'undetermined'

        if (temp_subtype, temp_lineage) in self.patterns:  #look for pattern from GISAID fasta file
            match = self.patterns[(temp_subtype, temp_lineage)]
            v['vtype'] = match[0].lower()
            v['subtype'] = match[1].lower()
            v['lineage'] = match[2].lower()
        elif vtype is not None and subtype is not None and lineage is not None: # defined from input for all viruses in input file
            v['vtype'] = vtype.lower()
            v['subtype'] = subtype.lower()
            v['lineage'] = lineage.lower()
        else:
            if '/' in v['strain']:
                #determine virus type from strain name
                fields = v['strain'].split('/')
                if re.match(r'^[abcd]$',fields[0].strip()):
                    v['vtype'] = fields[0].strip().lower()
                elif vtype is not None:
                    v['vtype'] = vtype.lower()
                else:
                    print("Couldn't parse virus type from strain name: " + v['strain'], (v['subtype'],v['lineage']))
            try:  # attempt to align with sequence from outgroup to determine subtype and lineage
                if v['strain'] not in self.db_strains:
                    scores = []
                    for olineage, oseq in self.outgroups.items():
                        SeqIO.write([oseq, seq], "temp_in.fasta", "fasta")
                        os.system("mafft --auto temp_in.fasta > temp_out.fasta 2>tmp")
                        tmp_aln = np.array(AlignIO.read('temp_out.fasta', 'fasta'))
                        scores.append((olineage, (tmp_aln[0]==tmp_aln[1]).sum()))
                    scores.sort(key = lambda x:x[1], reverse=True)
                    if scores[0][1]>0.85*len(seq):
                        print(v['strain'], (temp_subtype, temp_lineage), len(seq), "\n\t lineage based on similarity:",scores[0][0],"\n\t",scores)
                        match = self.outgroup_patterns[scores[0][0]]
                        v['vtype'] = match[0].lower()
                        v['subtype'] = match[1].lower()
                        v['lineage'] = match[2].lower()
            except:
                print("Couldn't parse virus subtype and lineage from aligning sequence: " + v['strain'], (temp_subtype, temp_lineage))

if __name__=="__main__":
    args = parser.parse_args()
    fasta_fields = {0: 'strain', 1: 'accession', 2: 'subtype', 4:'lineage', 5: 'date', 8: 'locus'}
    #              >B/Sao Paulo/61679/2014 | EPI_ISL_164700 | B / H0N0 | Original Sample | Yamagata | 2014-07-11 | Instituto Adolfo Lutz | Instituto Adolfo Lutz | HA
    #              >A/Alaska/81/2015 | EPI_ISL_197800 | A / H3N2 | Original |  | 2015-07-02 | Alaska State Virology Lab | Centers for Disease Control and Prevention | HA
    #              >A/Peru/06/2015 | EPI_ISL_209055 | A / H1N1 | C2/C1, MDCK1 | pdm09 | 2015-04-29 | Centers for Disease Control and Prevention | WHO Collaborating Centre for Reference and Research on Influenza | HA
    setattr(args, 'fasta_fields', fasta_fields)
    if args.path is None:
        args.path = "vdb/data/" + args.virus + "/"
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = Flu_vdb_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)