'''
Use this script to reconstruct titer database from archived files containing titer tables
Assumed directory structure is:
.
+-- data
|   +-- cdc
|       +-- subtype1
|           +-- flat_file_1.tsv
|           +-- flat_file_2.tsv
|           +-- ...
|       +-- subtype2
|       +-- ...
|   +-- elife
|       +-- subtype1
|           +-- flat_file_1.tsv
|           +-- flat_file_2.tsv
|           +-- ...
|       +-- subtype2
|       +-- ...
|   +-- nimr
|       +-- subtype1
|           +-- tabular_file_1.csv
|           +-- tabular_file_2.csv
|           +-- ...
|       +-- subtype2
|       +-- ...
+-- tdb
|   +-- upload_all.py
'''
import argparse
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument('-db', '--database', default='test_tdb', help="database to upload to")
parser.add_argument('--subtypes', nargs='+', type = str,  help ="flu subtypes to include, options are: h3n2, h1n1pdm, vic, yam")
parser.add_argument('--sources', nargs='+', type = str,  help ="data sources to include, options are: elife, nimr, cdc")
parser.add_argument('--nimr_path', default='data/nimr/', help="directory containing NIMR titers")
parser.add_argument('--cdc_path', default='data/cdc/', help="directory containing CDC titers")
parser.add_argument('--elife_path', default='data/elife/', help="directory containing eLife titers")

def upload_nimr(database, nimr_path, subtype):
    '''
    Makes calls to tdb/upload.py for every file in an NIMR titer directory (default data/nimr/).
    All files in the directory should be tabular titer tables in CSV format for only one subtype
    of virus (H3N2, H1N1pdm, etc.).
    '''
    path = nimr_path + subtype + "/"
    print "Uploading NIMR data for subtype", subtype, "contained in directory", path + "."

    for fname in os.listdir(path):
        if fname[0] != '.':
            fpath = path + fname
            fstem = fname[:-4]
            print "Uploading " + fname
            command = "python tdb/nimr_upload.py -db " + database + " --subtype " + subtype + " --ftype tables --path " + path + " --fstem " + fstem
            print "Running with: " + command
            subprocess.call(command, shell=True)
            print "Done with", fname + "."

    print "Done uploading NIMR documents."

def upload_cdc(database, cdc_path, subtype):
    '''
    Makes calls to tdb/upload.py for every flat file in an CDC titer directory (default data/cdc/).
    All files in the directory should be flat titer files in TSV format for only one subtype
    of virus (H3N2, H1N1pdm, etc.).
    '''
    path = cdc_path + subtype + "/"
    print "Uploading CDC data for subtype", subtype, "contained in directory", path + "."

    for fname in os.listdir(path):
        if fname[0] != '.':
            fpath = path + fname
            fstem = fname[:-4]
            print "Uploading " + fname
            command = "python tdb/cdc_upload.py -db " + database + " --subtype " + subtype + " --path " + path + " --fstem " + fstem
            print "Running with: " + command
            subprocess.call(command, shell=True)
            print "Done with", fname + "."

    print "Done uploading CDC documents."

def upload_elife(database, elife_path, subtype):
    '''
    Makes calls to tdb/upload.py for every flat file in an eLife titer directory (default data/elife/).
    All files in the directory should be flat titer files in TSV format for only one subtype
    of virus (H3N2, H1N1pdm, etc.).
    '''
    path = elife_path + subtype + "/"
    print "Uploading eLife data for subtype", subtype, "contained in directory", path + "."

    for fname in os.listdir(path):
        if fname[0] != '.':
            fpath = path + fname
            fstem = fname[:-4]
            print "Uploading " + fname
            command = "python tdb/elife_upload.py -db " + database + " --subtype " + subtype + " --path " + path + " --fstem " + fstem
            print "Running with: " + command
            subprocess.call(command, shell=True)
            print "Done with", fname + "."

    print "Done uploading stored eLife documents."

if __name__=="__main__":
    params = parser.parse_args()

    print "Beginning construction of", params.database + "."

    if params.subtypes is None:
        params.subtypes = ['h3n2', 'h1n1pdm', 'vic', 'yam']

    if params.sources is None:
        params.sources = ['elife', 'nimr', 'cdc']

    for source in params.sources:
        if source == "cdc":
            upload_cdc(params.database, params.cdc_path)
        if source == "elife":
            for subtype in params.subtypes:
                upload_elife(params.database, params.elife_path, subtype)
        if source == "nimr":
            for subtype in params.subtypes:
                upload_nimr(params.database, params.nimr_path, subtype)

    print "Done with all uploads."
