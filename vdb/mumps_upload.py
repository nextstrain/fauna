import os, re, time, datetime, csv, sys
from rethinkdb import r
from Bio import SeqIO
from upload import upload
from upload import get_parser

class mumps_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)
        self.strain_fix_fname = "source-data/mumps_strain_name_fix.tsv"
        self.location_fix_fname = "source-data/mumps_location_fix.tsv"
        self.date_fix_fname = "source-data/mumps_date_fix.tsv"

    def fix_name(self, name):
        original_name = name
        name = self.replace_strain_name(original_name, self.fix_whole_name)
        name = name.replace('MuV/', '').replace('MuVi/', '').replace('MuVs/','')
        name = name.replace('&','_').replace('?','_')
        name = re.sub(r'[_ ]?\[([A-Z])\]$', r'/\1', name)
        name = re.sub(r'\(([A-Z])\)$', r'/\1', name)
        name = re.sub(r'_([A-Z])_$', r'/\1', name)
        name = re.sub(r'[ ;]', r'_', name)
        name = re.sub(r'//', r'/', name)
        name = self.replace_strain_name(name, self.fix_whole_name)
        return name, original_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])

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
            self.define_MuV_genotype(doc)
            self.format_place(doc)
            self.format_region(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)

    def define_MuV_genotype(self, v):
    	MuV_genotypes_list = ["A","B","C","D","F","G","H","I","J","K","L","N"]
    	strain_name = v['strain']
    	if 'MuV_genotype' in v:
    		genotype = v['MuV_genotype']
    	else:
    		genotype = ""
    	if genotype == "" and strain_name != "" and strain_name.split("/")[-1] in MuV_genotypes_list:
    		MuV_genotype = strain_name.split("/")[-1]
    		v['MuV_genotype'] = MuV_genotype
    	return(v)

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    virus_fasta_fields = {1:'strain', 2:'collection_date', 3: 'host', 4:'country', 5:'division', 6: 'MuV_genotype'}
    sequence_fasta_fields = {0:'accession', 1:'strain'}
    # 0                          1                                2          3     4   5             6
    #>Massachusetts_outbreak_123|MuVs/Massachusetts.USA/50.16/1/G|2016-12-14|Human|USA|Massachusetts|G
    # 0       1                                  2         3     4      5                6
    #>BCCDC90|MuVs/BritishColumbia.CAN/34.16/1/G|2016-8-19|human|canada|british columbia|G
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = mumps_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
