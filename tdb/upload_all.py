#
# To use this upload script, fill in the paths containing NIMR, CDC, eLife, and stored
# titer tables under the appropriate variables in main marked with "TODO:"
#

import subprocess
import os
import shutil

def upload_nimr(database, directory):
    # Upload NIMR tables (September Report)
    print "Beginning upload of NIMR documents to", database + "."
    subtypes = ["h3n2", "h1n1pdm", "vic", "yam"]
    for subtype in subtypes:
        path = "data/" + subtype + "/"
        print "Uploading NIMR reports for subtype", subtype, "contained in directory", path + "."

        for fname in os.listdir(path):
            if fname [0] == "N":
                fpath = directory + path + fname
                fstem = fname[:-4]
                command = "python tdb/nimr_upload.py -db " + database + " --subtype " + subtype + " --ftype tables --path " + path + " --fstem " + fstem
                subprocess.call(command, shell=True)
                print "Done with", fname + "."

    print "Done uploading NIMR documents."

def upload_cdc(database, directory):
    # Upload CDC flat files
    print "Beginning upload of CDC documents to", database + "."
    subtypes = ["h3n2"]
    for subtype in subtypes:
        path = directory
        print "Uploading CDC reports for subtype", subtype, "contained in directory", path + "."

        for fname in os.listdir(path):
            if fname [0:3] == "test":
                fpath = path + fname
                fstem = fname[:-4]
                command = "python tdb/cdc_upload.py -db " + database + " --subtype " + subtype + " --path " + path + "tmp/ --fstem " + fstem
                subprocess.call(command, shell=True)
                print "Done with", fname + "."

    print "Done uploading CDC documents."

def upload_stored(database, directory):
    # Upload stored NIMR files
    print "Beginning upload of stored NIMR documents to", database + "."
    subtypes = ["h3n2", "h1n1pdm", "vic", "yam"]
    for subtype in subtypes:
        path = directory + subtype + "/"
        print "Uploading CDC reports for subtype", subtype, "contained in directory", path + "."

        for fname in os.listdir(path):
            if fname [0] == "N":
                fpath = path + fname
                fstem = fname[:-4]
                command = "python tdb/nimr_upload.py -db " + database + " --subtype " + subtype + " --ftype tables --path " + path + "tmp/"
                subprocess.call(command, shell=True)
                print "Done with", fname + "."

    print "Done uploading stored NIMR documents."

def upload_elife(database, directory):
    # Upload stored elife files
    print "Beginning upload of stored NIMR documents to", database + "."
    subtypes = ["h3n2", "vic", "yam"]
    for subtype in subtypes:
        path = directory + subtype + "/"
        print "Uploading CDC reports for subtype", subtype, "contained in directory", path + "."

        for fname in os.listdir(path):
            if fname [0] == "b":
                fpath = path + fname
                fstem = fname[:-4]
                command = "python tdb/elife_upload.py -db " + database + " --subtype " + subtype + " --path " + path + "tmp/ --fstem " + fstem
                subprocess.call(command, shell=True)
                print "Done with", fname + "."

    print "Done uploading stored eLife documents."



if __name__=="__main__":
    db = #TODO: Fill this in.
    nimr_dir = #TODO: Fill this in.
    cdc_dir = #TODO: Fill this in.
    stored_dir = #TODO: Fill this in.
    elife_dir = #TODO: Fill this in.
    print "Beginning construction of", db + "."
    upload_nimr(db, nimr_dir)
    upload_cdc(db, cdc_dir)
    upload_stored(db, stored_dir)
    upload_elife(db, elife_dir)
    print "Done with all uploads."
