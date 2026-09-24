import os, re, time, datetime, csv, sys
import rethinkdb as r
from Bio import SeqIO
from upload import upload
from upload import get_parser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def add_subtype_to_virus_fastas():
    #Append virus subtype to fasta fields, to ultimately create column in fauna
    input_fasta = "data/human_cov_genome.fasta"
    output_fasta = str(input_fasta.replace('.fasta',''))+"_annotated.fasta"
    #subtypes to look for
    cov_subtypes = ['OC43', 'HKU1', 'NL63', '229E']
    cov_types = {'OC43':'beta', 'HKU1':'beta', 'NL63':'alpha', '229E':'alpha'}

    sequences = []

    with open(input_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            #Fasta fields format: 'gb-id|strain|segment|date|host|country|subtype|virus species'
            new_record_list = record.description.split('|')

            #Annotate subtypes
            new_record_list = new_record_list+['None', 'None', 'full']
            for subtype in cov_subtypes:
                if subtype in record.description:
                    new_record_list[-3] = subtype.lower()
                    new_record_list[-2] = cov_types[subtype]
            #New fasta fields format: 'gb-id|strain|segment|date|host|country|subtype|virus species|subtype|type|sequence_locus'
            new_record_description = '|'.join(new_record_list)


            sequences.append(SeqRecord(record.seq, id=new_record_description, description=new_record_description))

    SeqIO.write(sequences, output_fasta, "fasta")

def add_subtype_to_gene_fastas():
    #Add virus subtype info to Spike and HE sequence files as well
    #Combine all Spike sequences into one file with virus type annotated
    genes = ['spike', 'he']
    cov_subtypes = ['OC43', 'HKU1', 'NL63', '229E']
    cov_types = {'NL63':'alpha', '229E':'alpha', 'OC43':'beta', 'HKU1':'beta'}


    for gene in genes:
        sequences = []
        for virus_subtype in cov_subtypes:
            input_fasta = 'data/'+virus_subtype+'_'+gene+'_genomealign.fasta'
            try:
                with open(input_fasta, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        #Fasta fields format: 'gb-id|strain|segment|date|host|country|subtype|virus species'
                        new_record_list = record.description.split('|')

                        #Change accesion number and strain name to indicate gene, for fauna
                        new_record_list[0] = new_record_list[0]+'_'+str(gene)
                        new_record_list[1] = new_record_list[1]+'_'+str(gene)

                        #Annotate subtypes
                        virus_type = cov_types[virus_subtype]
                        new_record_list = new_record_list+[virus_subtype.lower(), virus_type, gene]

                        #New fasta fields format: 'gb-id|strain|segment|date|host|country|subtype|virus species|subtype|type|sequence_locus'
                        new_record_description = '|'.join(new_record_list)

                        sequences.append(SeqRecord(record.seq, id=new_record_description,
                                                   description=new_record_description))
            except:
                pass

        output_fasta = 'data/human_cov_'+str(gene)+"_annotated.fasta"

        SeqIO.write(sequences, output_fasta, "fasta")


def combine_all_loci_fastas(vipr_upload, gene_upload):
#Combine full_seq, spike, he into one fasta file for upload
    if vipr_upload==False and gene_upload==False:
        loci = ['genome','spike','he']
    elif vipr_upload==False and gene_upload==True:
        loci = ['spike','he']
    elif vipr_upload==True and gene_upload==False:
        loci = ['genome']

    sequences = []
    for locus in loci:
        input_fasta = 'data/human_cov_'+locus+'_annotated.fasta'
        with open(input_fasta, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append(record)

    output_fasta = 'data/human_cov_annotated.fasta'

    SeqIO.write(sequences, output_fasta, "fasta")


class seasonal_corona_upload(upload):
    def __init__(self, **kwargs):
        upload.__init__(self, **kwargs)


    def fix_name(self, name):
        original_name = name
        try:
            name = 'V' + str(int(name))
        except:
            pass
        return name, original_name

    def fix_casing(self, document):
        for field in ['host']:
            if field in document and document[field] is not None:
                document[field] = self.camelcase_to_snakecase(document[field])


if __name__=="__main__":
    parser = get_parser()
    parser.add_argument('--vipr_upload', action='store_true', help='specify whether to upload sequences downloaded from Vipr')
    parser.add_argument('--gene_upload', action='store_true', help='specify whether to upload sequences Spike and HE sequences extracted from alignments')
    args = parser.parse_args()

    if args.vipr_upload == True and args.gene_upload == False:
        add_subtype_to_virus_fastas()
    if args.gene_upload == True and args.vipr_upload == False:
        add_subtype_to_gene_fastas()
    if args.vipr_upload ==False and args.gene_upload ==False:
        add_subtype_to_virus_fastas()
        add_subtype_to_gene_fastas()

    combine_all_loci_fastas(vipr_upload = args.vipr_upload, gene_upload = args.gene_upload)

    #gb-id|strain|segment|date|host|country|subtype|virus species|subtype|type
    virus_fasta_fields = {1:'strain', 3:'collection_date', 4: 'host', 5:'country', 7:'virus_species', 8:'subtype', 9:'type', 10: 'sequence_locus'}
    sequence_fasta_fields = {0:'accession', 1:'strain', 8:'subtype', 9:'type', 10: 'sequence_locus'}
    setattr(args, 'virus_fasta_fields', virus_fasta_fields)
    setattr(args, 'sequence_fasta_fields', sequence_fasta_fields)
    connVDB = seasonal_corona_upload(**args.__dict__)
    connVDB.upload(**args.__dict__)
