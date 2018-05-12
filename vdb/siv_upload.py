import re
from upload import upload
from upload import get_parser

class siv_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)

    def fix_name(self, name):
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('//', '/').replace('.', '').replace(',', '')
        try:
            tmp_name = 'SIV' + str(int(tmp_name))
        except:
            pass
        return tmp_name

    def filter(self, documents, index, **kwargs):
        for doc in documents:
            if doc[index] == "":
                print(doc[index])
        result = filter(lambda doc: doc[index] != "" and doc[index] is not None,documents)
        return result

    def format(self, documents, exclude_virus_methods=False, **kwargs):
        '''
        format virus information in preparation to upload to database table
        '''
        self.define_regions("source-data/geo_regions.tsv")
        self.define_countries("source-data/geo_synonyms_siv.tsv")
        self.countries_siv = set()
        for doc in documents:
            if 'strain' in doc:
                doc['strain'] = self.fix_name(doc['strain'])
            self.format_date(doc)
            self.format_place(doc)
            self.format_country(doc)
            self.format_region(doc)
            self.rethink_io.check_optional_attributes(doc, [])
            self.fix_casing(doc)

    def format_country(self, v):
        '''
        Label viruses with country based on strain name
        '''
        if 'country' in v and v['country'] != '' and v['country'] is not None:
            original_country = v['country']
            result = self.determine_location(v['country'])
            if result is not None:
                v['country'], v['division'], v['location'] = result
            else:
                v['country'], v['division'], v['location'] = None, None, None
                print("couldn't parse country for ",  v['strain'], original_country)
        else:
            v['country'], v['division'], v['location'] = None, None, None
            #print(v['strain'], " country field isn't defined")

    def determine_location(self, name):
        '''
        Try to determine country, division and location information from name
        Return tuple of country, division, location if found, otherwise return None
        '''
        try:
            label = re.match(r'^([^/]+)', name).group(1).lower()						# check first for whole geo match
            if label in self.label_to_country:
                return (self.label_to_country[label], self.label_to_division[label], self.label_to_location[label])
            else:
                return None
        except:
            return None

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    virus_fasta_fields = {2:'host_species', 3:'sub_species', 4:'SIVxyz', 5:'strain', 8:'country', 9:'collection_date'}
    sequence_fasta_fields = {0:'accession', 1:'LanL_ID', 6:'sequence_length'}
    # 0        1      2  3          4     5         6  7
    #>>accession|lanl_ID|host_species|subspecies|SIVxyz|lanl name_isolate name|sequence_length|region|country|year
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = siv_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
