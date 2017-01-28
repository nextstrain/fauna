import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument('--db', '--database', default='test_tdb_2', help="database to upload to")
parser.add_argument('--subtype', default='h3n2', help="subtype to upload")
parser.add_argument('--nimr_path', default='data/nimr/', help="directory containing NIMR titers")
parser.add_argument('--cdc_path', default='data/cdc/', help="directory containing CDC titers")
parser.add_argument('--elife_path', default='data/elife/', help="directory containing eLife titers")

def upload_nimr(database, directory, nimr_path, subtype):
    # Upload NIMR tables
    print "Beginning upload of NIMR documents to", database + "."
    path = nimr_path + subtype + "/"
    print "Uploading NIMR reports for subtype", subtype, "contained in directory", path + "."

    for fname in os.listdir(path):
        fpath = path + fname
        fstem = fname[:-4]
        command = "python tdb/nimr_upload.py -db " + database + " --subtype " + subtype + " --ftype tables --path " + path + " --fstem " + fstem
        subprocess.call(command, shell=True)
        print "Done with", fname + "."

    print "Done uploading NIMR documents."

def upload_cdc(database, cdc_path, subtype):
    # Upload CDC flat files
    print "Beginning upload of CDC documents to", database + "."
    path = cdc_path + subtype + "/"
    print "Uploading CDC reports for subtype", subtype, "contained in directory", path + "."

    for fname in os.listdir(path):
        fpath = path + fname
        fstem = fname[:-4]
        command = "python tdb/cdc_upload.py -db " + database + " --subtype " + subtype + " --path " + path + " --fstem " + fstem
        subprocess.call(command, shell=True)
        print "Done with", fname + "."

    print "Done uploading CDC documents."

def upload_elife(database, elife_path, subtype):
    # Upload stored elife files
    print "Beginning upload of stored NIMR documents to", database + "."
    path = elife_path + subtype + "/"
    print "Uploading CDC reports for subtype", subtype, "contained in directory", path + "."

    for fname in os.listdir(path):
        fpath = path + fname
        fstem = fname[:-4]
        command = "python tdb/elife_upload.py -db " + database + " --subtype " + subtype + " --path " + path + " --fstem " + fstem
        subprocess.call(command, shell=True)
        print "Done with", fname + "."

    print "Done uploading stored eLife documents."

if __name__=="__main__":
    args = parser.parse_args()
    print "Beginning construction of", db + "."
    upload_nimr(db, args.nimr_dir, "h3n2")
    upload_nimr(db. args.nimr_dir, "h1n1pdm")
    upload_nimr(db. args.nimr_dir, "vic")
    upload_nimr(db. args.nimr_dir, "yam")
    upload_cdc(db, args.cdc_dir, "h3n2")
    upload_elife(db, args.elife_dir, "h3n2")
    upload_elife(db, args.elife_dir, "vic")
    upload_elife(db, args.elife_dir, "yam")
    print "Done with all uploads."
