import os, json, datetime, sys
import rethinkdb as r
sys.path.append('')  # need to import from base
from base.rethink_io import rethink_io
from vdb.download import download as vdb_download

def get_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-db', '--database', default='vdb', help="database to download from")
    parser.add_argument('--rethink_host', default=None, help="rethink host url")
    parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
    parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
    parser.add_argument('-v', '--virus', help="virus name")
    parser.add_argument('--subtype', default='h3n2', help="subtype to be include in download")
    parser.add_argument('--ftype', default='fasta', help="output file format, default \"fasta\", other options are \"json\" and \"tsv\"")
    parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
    parser.add_argument('--path', default='data', help="path to dump output files to")
    parser.add_argument('--fasta_fields', default=['strain', 'virus', 'accession', 'date', 'region', 'country', 'division', 'location', 'source', 'locus', 'authors'], help="fasta fields for output fasta")

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
        self.database = database.lower()
        self.measurements = []

    def connect_rethink(self, **kwargs):
        if self.database not in ['tdb', 'test_tdb', 'test', 'test_tdb_2']:
            raise Exception("Cant download to this database: " + self.database)
        self.rethink_io = rethink_io()
        self.rethink_host, self.auth_key = self.rethink_io.assign_rethink(**kwargs)
        self.rethink_io.connect_rethink(self.database, self.rethink_host, self.auth_key)
        self.rethink_io.check_table_exists(self.database, self.virus)

    def count_documents(self):
        '''
        return integer count of number of documents in table
        '''
        return r.db(self.database).table(self.virus).count().run()

    def download(self, subtype, output=True, count=True, **kwargs):
        '''
        download documents from table
        '''
        import time
        start_time = time.time()
        self.connect_rethink(**kwargs)
        self.vdb_download = vdb_download(database=self.database, virus=self.virus)
        select, present, interval, = self.vdb_download.parse_subset_arguments(**kwargs)
        select.append(('subtype', [subtype]))
        sequence_count = r.table(self.virus).count().run()
        print(sequence_count, "measurements in table:", self.virus)
        print("Downloading titer measurements from the table: " + self.virus)
        measurements = self.rethinkdb_download(self.virus, presents=present, selections=select, intervals=interval, **kwargs)
        print("Downloaded " + str(len(measurements)) + " measurements")
        if output:
            self.output(measurements, subtype=subtype, **kwargs)
        if count:
            self.write_count(measurements, subtype=subtype, **kwargs)
        print("--- %s minutes to download ---" % ((time.time() - start_time)/60))

    def rethinkdb_download(self, table, **kwargs):
        '''
        Default command merges documents from the sequence table and virus table
        Chain rethinkdb filter and has_fields commands to the default command
        Return documents from the database that are left after filtering
        '''
        # take each sequence and merge with corresponding virus document
        command = r.table(table)
        command = self.vdb_download.add_present_command(command, **kwargs)
        command = self.vdb_download.add_selections_command(command, **kwargs)
        command = self.vdb_download.add_intervals_command(command, **kwargs)
        sequences = list(command.run())
        return list(sequences)

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
            print("Wrote to " + fname)

    def write_text(self, measurements, fname, text_fields=['virus_strain', 'serum_strain', 'serum_id', 'source', 'titer']):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for meas in measurements:
                for field in text_fields:
                    if field in meas and meas[field] is not None:
                        handle.write(str(meas[field]) + '\t')
                    else:
                        handle.write('?' + '\t')
                handle.write('\n')
            handle.close()
            print("Wrote to " + fname)

    def output(self, measurements, path, fstem, ftype, subtype, **kwargs):
        fname = path + '/' + fstem + '.' + ftype
        if ftype == 'json':
            self.write_json(measurements,fname)
        elif ftype == 'tsv':
            self.write_text(measurements, fname)
        elif ftype == 'augur':
            fname = path + '/' + subtype + "_hi_titers.tsv"
            self.write_text(measurements, fname)
        else:
            raise Exception("Can't output to that file type, only json or text allowed")

    def write_count(self, measurements, path, subtype, **kwargs):
        fname = path + '/' + subtype + '_hi_strains.tsv'
        print("Counting HI_strains and printing to " + fname)
        HI_titer_count = self.count(measurements)
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            for strain, count in HI_titer_count.items():
                handle.write(strain + "\t" + str(count))
                handle.write("\n")
            handle.close()
            print("Wrote to " + fname)

    def count(self, measurements):
        HI_titer_count = {}
        for meas in measurements:
            if meas['serum_strain'] in HI_titer_count.keys():
                HI_titer_count[meas['serum_strain']] = HI_titer_count[meas['serum_strain']] + 1
            else:
                HI_titer_count[meas['serum_strain']] = 1
            if meas['virus_strain'] in HI_titer_count.keys():
                HI_titer_count[meas['virus_strain']] = HI_titer_count[meas['virus_strain']] + 1
            else:
                HI_titer_count[meas['virus_strain']] = 1
        return HI_titer_count

if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()
    current_date = str(datetime.datetime.strftime(datetime.datetime.now(),'%Y_%m_%d'))
    if args.fstem is None:
        args.fstem = args.subtype + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    connTDB = download(**args.__dict__)
    connTDB.download(**args.__dict__)
