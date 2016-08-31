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
    parser.add_argument('--interval', nargs='+', type=str, default=[], help="Select interval of values for fields \'--interval field1:value1,value2 field2:value1,value2\'")
    parser.add_argument('--years_back', type=str, default=None, help='number of past years to sample sequences from \'--years_back field:value\'')
    parser.add_argument('--relaxed_interval', default=False, action="store_true", help="Relaxed comparison to date interval, 2016-XX-XX in 2016-01-01 - 2016-03-01")
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

    def connect_rethink(self, **kwargs):
        if self.database not in ['vdb', 'test_vdb', 'test']:
            raise Exception("Can't download from this database: " + self.database)
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

        import time
        start_time = time.time()
        self.connect_rethink(**kwargs)
        select, present, interval= self.parse_subset_arguments(**kwargs)
        sequence_count = r.table(self.sequences_table).count().run()
        print(sequence_count, "sequences in table:", self.sequences_table)
        virus_count = r.table(self.viruses_table).count().run()
        print(virus_count, "viruses in table:", self.viruses_table)
        print("Downloading documents from the sequence table: " + self.sequences_table, " and virus table: ", self.viruses_table)
        sequences = self.rethinkdb_download(self.sequences_table, self.viruses_table, presents=present, selections=select, intervals=interval, **kwargs)
        print("Downloaded " + str(len(sequences)) + " sequences")
        sequences = self.resolve_duplicates(sequences, **kwargs)
        if output:
            self.output(sequences, **kwargs)
        print("--- %s minutes to download ---" % ((time.time() - start_time)/60))

    def parse_subset_arguments(self, select=[], present=[], interval=[], years_back=None, **kwargs):
        '''
        Parse arguments needed for subsetting on the server side
        Return a tuple of the arguments
        '''
        selections = self.parse_select_argument(select)
        intervals = self.parse_select_argument(interval)
        if years_back is not None and len(intervals)==0:
            intervals = self.parse_years_back_argument(years_back)
        return selections, present, intervals,

    def parse_select_argument(self, groupings=[]):
        '''
        :arg groupings like country:brazil,argentina
        parse the 'select' parameter to determine which field name to filter and for what values
        :return: [(grouping name, [group values])] ex. [(country, [brazil, argentina)]
        '''
        selections = []
        if groupings is not None and len(groupings)>0:
            for group in groupings:
                result = group.split(':')
                selections.append((result[0].lower(), result[1].lower().split(',')))
        return selections

    def parse_years_back_argument(self, argument):
        field = argument.split(':')[0].lower()
        years = int(argument.split(':')[1])
        date_now = datetime.datetime.utcnow()
        older_date = str(datetime.datetime.strftime(date_now.replace(year=date_now.year-years),'%Y-%m-%d'))
        newer_date = str(datetime.datetime.strftime(date_now,'%Y-%m-%d'))
        return [(field, [older_date, newer_date])]

    def rethinkdb_download(self, sequence_table, virus_table, index='strain', **kwargs):
        '''
        Default command merges documents from the sequence table and virus table
        Chain rethinkdb filter and has_fields commands to the default command
        Return documents from the database that are left after filtering
        '''
        # take each sequence and merge with corresponding virus document
        command = r.table(sequence_table).merge(lambda sequence: r.table(virus_table).get(sequence[index]))
        command = self.add_present_command(command, **kwargs)
        command = self.add_selections_command(command, **kwargs)
        command = self.add_intervals_command(command, **kwargs)
        command = self.add_public_command(command, **kwargs)
        sequences = list(command.run())
        return list(sequences)

    def add_present_command(self, command, presents=[], **kwargs):
        '''
        Add present fields check to command
        '''
        if len(presents)>0:
            print("Only downloading documents with fields: " + str(presents))
            command = command.has_fields(r.args(presents))
        return command

    def add_selections_command(self, command, selections=[], **kwargs):
        '''
        Add selections filter to command
        '''
        if len(selections)>0:
            for sel in selections:
                field = sel[0]
                values = sel[1]
                print("Only downloading documents with field \'" + field + "\' equal to one of " + str(values))
                command = command.filter(lambda doc: r.expr(values).contains(doc[field]))
        return command

    def add_intervals_command(self, command, intervals=[], relaxed_interval=False, **kwargs):
        '''
        Add intervals filter to command
        '''
        if len(intervals) > 0:
            for sel in intervals:
                field = sel[0]
                values = sel[1]
                older_date, newer_date = self.check_date_format(values[0], values[1])
                command = command.filter(lambda doc: rethinkdb_date_greater(doc[field].split('-'), older_date.split('-'), relaxed_interval))
                command = command.filter(lambda doc: rethinkdb_date_greater(newer_date.split('-'), doc[field].split('-'), relaxed_interval))
                print('Only downloading documents in the interval specified (' + ' - '.join([older_date, newer_date]) + ') for field \'' + field + '\'')
                if relaxed_interval:
                    print("Using relaxed interval comparison")
        return command

    def add_public_command(self, command, public_only, **kwargs):
        '''
        Add public filter to command
        '''
        if public_only:
            print("Only downloading public sequences")
            command = command.filter({'public': True})
        return command

    def check_date_format(self, older_date, newer_date):
        one_sided_symbols = ['', 'XXXX-XX-XX']
        if newer_date in one_sided_symbols:
            newer_date = str(datetime.datetime.strftime(datetime.datetime.utcnow(),'%Y-%m-%d'))
        if older_date in one_sided_symbols:
            older_date = '0000-00-00'
        if older_date > newer_date:
            raise Exception("Date interval must list the earlier date first")
        if not re.match(r'\d\d\d\d-(\d\d)-(\d\d)$', older_date) or not re.match(r'\d\d\d\d-(\d\d)-(\d\d)$', newer_date):
            raise Exception("Date interval must be in YYYY-MM-DD format with all values defined", older_date, newer_date)
        return(older_date.upper(), newer_date.upper())

    def resolve_duplicates(self, sequences, pick_longest=True, **kwargs):
        strain_locus_to_doc = {doc['strain']+doc['locus']: doc for doc in sequences}
        if pick_longest:
            print("Resolving duplicate strains and locus by picking the longest sequence")
            for doc in sequences:
                if doc['strain']+doc['locus'] in strain_locus_to_doc:
                    if self.longer_sequence(doc['sequence'], strain_locus_to_doc[doc['strain']+doc['locus']]):
                        strain_locus_to_doc[doc['strain']+doc['locus']] = doc
                else:
                    strain_locus_to_doc[doc['strain']] = doc
        return list(strain_locus_to_doc.values())

    def longer_sequence(self, long_seq, short_seq):
        '''
        :return: true if long_seq is longer than short_seq
        '''
        return len(long_seq) > len(short_seq)

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

    def write_fasta(self, viruses, fname, sep='|', fasta_fields=['strain', 'virus', 'accession'], **kwargs):
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

    def write_tsv(self, viruses, fname, sep='\t', fasta_fields=['strain', 'virus', 'accession'], **kwargs):
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

    def output(self, documents, path, fstem, ftype, **kwargs):
        fname = path + '/' + fstem + '.' + ftype
        print("Outputing", len(documents), "documents to ", fname)
        if ftype == 'json':
            self.write_json(documents,fname)
        elif ftype == 'fasta':
            self.write_fasta(documents, fname, **kwargs)
        elif ftype == 'tsv':
            self.write_tsv(documents, fname, **kwargs)
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

def rethinkdb_date_greater(greater_date, comparison_date, relaxed_interval):
    return r.branch(r.lt(greater_date[0], comparison_date[0]),
        False,
        r.eq(greater_date[0], comparison_date[0]),
        r.branch(r.eq(greater_date[1], 'XX').or_(r.eq(comparison_date[1],'XX')),
            relaxed_interval,
            r.lt(greater_date[1], comparison_date[1]),
            False,
            r.eq(greater_date[1], comparison_date[1]),
            r.branch(r.eq(greater_date[2], 'XX').or_(r.eq(comparison_date[2],'XX')),
                relaxed_interval,
                r.lt(greater_date[2], comparison_date[2]),
                False,
                True),
            True),
        True)