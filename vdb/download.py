import os, json, datetime, sys, re
import rethinkdb as r
from Bio import SeqIO
import numpy as np
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io

def get_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-db', '--database', default='vdb', help="database to download from")
    parser.add_argument('--rethink_host', default=None, help="rethink host url")
    parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
    parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
    parser.add_argument('-v', '--virus', help="virus name")
    parser.add_argument('--ftype', default='fasta', help="output file format, default \"fasta\", other options are \"json\" and \"tsv\"")
    parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
    parser.add_argument('--path', default='data', help="path to dump output files to")
    parser.add_argument('--fasta_fields', default=['strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'location', 'source', 'locus', 'authors'], help="fasta fields for output fasta")

    parser.add_argument('--public_only', default=False, action="store_true", help="include to subset public sequences")
    parser.add_argument('--select', nargs='+', type=str, default=[], help="Select specific fields ie \'--select field1:value1 field2:value1,value2\'")
    parser.add_argument('--present', nargs='+', type=str, default=[], help="Select specific fields to be non-null ie \'--present field1 field2\'")
    parser.add_argument('--interval', nargs='+', type=str, default=None, help="Select interval of values for fields \'--interval field1:value1,value2 field2:value1,value2\'")

    parser.add_argument('--pick_longest', default=False, action="store_true",  help ="For duplicate strains, only includes the longest sequence for each locus")
    parser.add_argument('--pick_recent', default=False, action="store_true",  help ="For duplicate strains, only includes the most recent sequence for each locus")
    return parser

class download(object):
    def __init__(self, database, virus, **kwargs):
        '''
        parser for virus, fasta fields, output file names, output file format path, interval
        '''
        self.virus = virus.lower()
        self.viruses_table = virus + "_viruses"
        self.sequences_table = virus + "_sequences"
        self.database = database.lower()
        if self.database not in ['vdb', 'test_vdb']:
            raise Exception("Cant download from this database: " + self.database)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.viruses_table)
        self.rethink_io.check_table_exists(self.database, self.sequences_table)

    def count_documents(self, table):
        '''
        return integer count of number of documents in table
        '''
        return r.db(self.database).table(table).count().run()

    def download(self, output=True, **kwargs):
        '''
        download documents from table
        '''
        print("Downloading all viruses from the table: " + self.viruses_table)
        viruses = list(r.table(self.viruses_table).run())
        print("Downloading all viruses from the table: " + self.sequences_table)
        sequences = list(r.table(self.sequences_table).run())
        sequences = self.link_viruses_to_sequences(viruses, sequences)
        sequences = self.subset(sequences, **kwargs)
        self.resolve_duplicates(sequences, **kwargs)
        if output:
            self.output(sequences, **kwargs)

    def resolve_duplicates(self, sequences, pick_longest, pick_recent, **kwargs):
        strain_locus_to_doc = {doc['strain']+doc['locus']: doc for doc in sequences}
        if pick_longest or pick_recent:
            for doc in sequences:
                if doc['strain']+doc['locus'] in strain_locus_to_doc:
                    if pick_longest and self.longer_sequence(doc['sequence'], strain_locus_to_doc[doc['strain']+doc['locus']]):
                        strain_locus_to_doc[doc['strain']] = doc
                    elif pick_recent and self.most_recent_sequence(doc['date'], strain_locus_to_doc[doc['strain']+doc['locus']]['date']):
                        strain_locus_to_doc[doc['strain']] = doc
                else:
                    strain_locus_to_doc[doc['strain']] = doc

    def subset(self, sequences, public_only=False, select=[], present=[], interval=[], **kwargs):
        selections = self.parse_select_argument(select)
        if public_only:
            selections.append(('public', True))
        intervals = self.parse_select_argument(interval)

        print("Documents from " + self.database + '.' + self.sequences_table + " before subsetting: " + str(len(sequences)))
        sequences = self.subsetting(sequences, selections, present, intervals, **kwargs)
        print("Documents from " + self.database + '.' + self.sequences_table + " after subsetting: " + str(len(sequences)))
        return sequences

    def subsetting(self, cursor, selections=[], presents=[], intervals=[], **kwargs):
        '''
        filter through documents in vdb to return subsets of sequence
        '''
        # determine fields in all documents that can be filtered by
        list_of_fields = [set(doc.keys()) for doc in cursor]
        fields = set.intersection(*list_of_fields)
        if len(selections) > 0:
            for sel in selections:
                if sel[0] in fields:
                    cursor = filter(lambda doc: doc[sel[0]] in sel[1], cursor)
                    print('Removed documents that were not in values specified (' + ','.join(sel[1]) + ') for field \'' + sel[0] + '\', remaining documents: ' + str(len(cursor)))
                else:
                    print(sel[0] + " is not in all documents, can't subset by that field")
        if len(presents) > 0:
            for sel in presents:
                if sel in fields:
                    cursor = filter(lambda doc: doc[sel] is not None, cursor)
                    print('Removed documents that were null for field \'' + sel + '\', remaining documents: ' + str(len(cursor)))
                else:
                    print(sel + " is not in all documents, can't subset by that field")
        if len(intervals) > 0:
            for sel in intervals:
                if sel[0] == 'collection_date' or sel[0] == 'submission_date':
                    older_date, newer_date = self.check_date_format(sel[1][0], sel[1][1])
                    cursor = filter(lambda doc: self.in_date_interval(doc[sel[0]], older_date, newer_date), cursor)
                    print('Removed documents that were not in the interval specified (' + ' - '.join(sel[1]) + ') for field \'' + sel[0] + '\', remaining documents: ' + str(len(cursor)))
                else:
                    print(sel[0] + " is not in all documents, can't subset by that field")
        return cursor

    def parse_select_argument(self, grouping):
        '''
        parse the 'select' parameter to determine which field name to filter and for what values
        :return: (grouping name, group values))
        '''
        selections = []
        if grouping is not None:
            for group in grouping:
                result = group.split(':')
                selections.append((result[0].lower(), result[1].lower().split(',')))
        return selections

    def check_date_format(self, older_date, newer_date):
        if newer_date == 'now' or newer_date == 'present':
            newer_date = str(datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y-%m-%d'))
        if older_date > newer_date:
            raise Exception("Date interval must list the earlier date first")
        if not re.match(r'\d\d\d\d-(\d\d)-(\d\d)', older_date) or not re.match(r'\d\d\d\d-(\d\d)-(\d\d)', newer_date):
            raise Exception("Date interval must be in YYYY-MM-DD format with all values defined")
        return(older_date, newer_date)

    def in_date_interval(self, date, older_date, newer_date):
        '''
        :return: true if the date is in the interval older_date - newer_date, otherwise False
        '''
        return self.date_greater(date.split('-'), older_date.split('-')) and self.date_greater(newer_date.split('-'), date.split('-'))

    def date_greater(self, greater_date, comparison_date):
        '''
        :return: true if greater_date > comparison_date
        '''
        # compare year
        if greater_date[0] < comparison_date[0]:
            return False
        elif greater_date[0] == comparison_date[0]:
            # compare month
            if greater_date != 'XX' and comparison_date != 'XX' and greater_date[1] < comparison_date[1]:
                return False
            elif greater_date[1] == comparison_date[1]:
                # compare day
                if greater_date != 'XX' and comparison_date != 'XX' and greater_date[2] < comparison_date[2]:
                    return False
        return True

    def link_viruses_to_sequences(self, viruses, sequences):
        '''
        copy the virus doc to the sequence doc
        '''
        linked_sequences = []
        strain_name_to_virus_doc = {virus['strain']: virus for virus in viruses}
        for sequence_doc in sequences:
            if sequence_doc['strain'] in strain_name_to_virus_doc:  # determine if sequence has a corresponding virus to link to
                virus_doc = strain_name_to_virus_doc[sequence_doc['strain']]
                sequence_doc.update(virus_doc)
                linked_sequences.append(sequence_doc)
        return linked_sequences

    def longer_sequence(self, long_seq, short_seq):
        '''
        :return: true if long_seq is longer than short_seq
        '''
        return long_seq > short_seq

    def most_recent_sequence(self, new_seq, old_seq):
        '''
        :return: true if new_seq is more recent than old_seq
        '''
        return new_seq > old_seq

    def write_json(self, data, fname, indent=1):
        '''
        writes as list of viruses (dictionaries)
        '''
        try:
            handle = open(fname, 'w')
        except:
            print("Couldn't open output file")
            print(fname)
            raise FileNotFoundError
        else:
            json.dump(data, handle, indent=indent)
            handle.close()

    def write_fasta(self, viruses, fname, sep='|', fasta_fields=['strain', 'virus', 'accession', 'date', 'region',
                                                                 'country', 'division', 'location', 'source', 'locus',
                                                                 'authors']):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for virus in viruses:
            	if 'sequence' in virus:
            	    if virus['sequence']:
                        fields = [str(virus[field]) if (field in virus and virus[field] is not None) else '?'
                                  for field in fasta_fields]
                        handle.write(">"+sep.join(fields)+'\n')
                        handle.write(virus['sequence'] + "\n")
            handle.close()

    def write_tsv(self, viruses, fname, sep='\t', fasta_fields=['strain', 'virus', 'accession', 'date', 'region',
                                                                 'country', 'division', 'location', 'source', 'locus',
                                                                 'authors']):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            handle.write(sep.join(fasta_fields)+'\n')
            for virus in viruses:
                fields = [str(virus[field]) if (field in virus and virus[field] is not None) else '?'
                          for field in fasta_fields]
                handle.write(sep.join(fields)+'\n')
            handle.close()


    def output(self, documents, path, fstem, ftype, fasta_fields, **kwargs):
        fname = path + '/' + fstem + '.' + ftype
        print("Outputing", len(documents), "documents to ", fname)
        if ftype == 'json':
            self.write_json(documents,fname)
        elif ftype == 'fasta':
            self.write_fasta(documents, fname, fasta_fields=fasta_fields+self.virus_specific_fasta_fields)
        elif ftype == 'tsv':
            self.write_tsv(documents, fname, fasta_fields=fasta_fields+self.virus_specific_fasta_fields)
        else:
            raise Exception("Can't output to that file type, only json, fasta or tsv allowed")
        print("Wrote to " + fname)

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connVDB = download(**args.__dict__)
    connVDB.download(**args.__dict__)