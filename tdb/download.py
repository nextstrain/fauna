import os, json, datetime, sys, time
from rethinkdb import r

# Enable import from modules in parent directory.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from base.rethink_io import rethink_io
from vdb.download import download as vdb_download

def get_parser():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-db', '--database', nargs='+', type=str, default=['tdb'], help="database(s) to download from")
    parser.add_argument('--rethink_host', default=None, help="rethink host url")
    parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
    parser.add_argument('--local', default=False, action="store_true",  help ="connect to local instance of rethinkdb database")
    parser.add_argument('-v', '--virus', default='flu', help="virus name")
    parser.add_argument('--subtype', help="subtype to be included in download")
    parser.add_argument('--ftype', default='tsv', help="output file format, default \"tsv\", options are \"json\" and \"tsv\"")
    parser.add_argument('--fstem', default=None, help="default output file name is \"VirusName_Year_Month_Date\"")
    parser.add_argument('--path', default='data', help="path to dump output files to")

    parser.add_argument('--select', nargs='+', type=str, default=[], help="Select specific fields ie \'--select field1:value1 field2:value1,value2\'")
    parser.add_argument('--present', nargs='+', type=str, default=[], help="Select specific fields to be non-null ie \'--present field1 field2\'")
    parser.add_argument('--interval', nargs='+', type=str, default=[], help="Select interval of values for fields \'--interval field1:value1,value2 field2:value1,value2\'")
    parser.add_argument('--years_back', type=str, default=None, help='number of past years to sample sequences from \'--years_back field:value\'')
    parser.add_argument('--relaxed_interval', default=False, action="store_true", help="Relaxed comparison to date interval, 2016-XX-XX in 2016-01-01 - 2016-03-01")
    return parser

class download(object):
    def __init__(self, database, virus, **kwargs):
        '''
        parser for virus, output file names, output file format path, interval
        '''
        self.virus = virus.lower()
        self.database = database.lower()
        self.measurements = []

    def connect_rethink(self, **kwargs):

        if self.database not in ['tdb', 'test_tdb', 'cdc_tdb', 'test_cdc', 'crick_tdb', 'test_crick', 'niid_tdb', 'test_niid', 'vidrl_tdb', 'test_vidrl']:
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

    def download(self, subtype=None, **kwargs):
        '''
        download documents from table
        '''
        start_time = time.time()
        self.connect_rethink(**kwargs)
        self.vdb_download = vdb_download(database=self.database, virus=self.virus)
        select, present, interval, = self.vdb_download.parse_subset_arguments(**kwargs)

        if subtype is not None:
            select.append(('subtype', [subtype]))

        measurement_count = r.table(self.virus).count().run()
        print(measurement_count, "measurements in table:", self.virus)
        print("Downloading titer measurements from the table: " + self.virus)
        measurements = self.rethinkdb_download(self.virus, presents=present, selections=select, intervals=interval, **kwargs)
        self.rename_strains_with_passage(measurements)
        print("Downloaded " + str(len(measurements)) + " measurements")
        print("--- %s minutes to download ---" % ((time.time() - start_time)/60))
        return measurements

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

    def rename_strains_with_passage(self, measurements):
        '''
        Append -egg to virus_strain and serum_strain for egg-passaged viruses
        '''
        for measurement in measurements:
            if measurement['virus_passage_category'] == "egg":
                measurement['virus_strain'] = measurement['virus_strain'] + "-egg"
            if measurement['serum_passage_category'] == "egg":
                measurement['serum_strain'] = measurement['serum_strain'] + "-egg"

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

    def write_text(self, measurements, fname, text_fields=['virus_strain', 'serum_strain', 'serum_id', 'source', 'titer', 'assay_type']):
        try:
            handle = open(fname, 'w')
        except IOError:
            pass
        else:
            # Write out field names as a header.
            handle.write("\t".join(text_fields) + "\n")

            # Write out field values for each record.
            for meas in measurements:
                for field in text_fields:
                    if field in meas and meas[field] is not None:
                        handle.write(str(meas[field]) + '\t')
                    else:
                        handle.write('?' + '\t')
                handle.write('\n')
            handle.close()
            print("Wrote to " + fname)

    def output(self, measurements, path, fstem, ftype, **kwargs):
        fname = path + '/' + fstem + '_titers' + '.' + ftype
        if ftype == 'json':
            self.write_json(measurements,fname)
        elif ftype == 'tsv':
            self.write_text(measurements, fname)
        else:
            raise Exception("Can't output to that file type, only json or text allowed")

    def write_count(self, measurements, path, fstem, **kwargs):
        fname = path + '/' + fstem + '_strains.tsv'
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
        if args.subtype is not None:
            args.fstem = args.subtype + '_' + current_date
        else:
            args.fstem = args.virus + '_' + current_date
    if not os.path.isdir(args.path):
        os.makedirs(args.path)
    measurements = []
    for database in args.database:
        print("Downloading from " + database + " database")
        dict_args = vars(args)
        dict_args["database"] = database
        print(dict_args)
        downloader = download(**dict_args)
        db_measurements = downloader.download(**dict_args)
        measurements.extend(db_measurements)
    downloader.output(measurements, **dict_args)
