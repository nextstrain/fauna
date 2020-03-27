from __future__ import print_function
import os, subprocess
import re

def fix_input_file_spaces():
    path = '../../fludata/VIDRL-Melbourne-WHO-CC/raw-data/'
    for ab in os.listdir(path):
        ab_path = '{}{}/'.format(path, ab)
        for subtype in os.listdir(ab_path):
            subtype_path = '{}{}/'.format(ab_path, subtype)
            for assay in os.listdir(subtype_path):
                complete_path = '{}{}/'.format(subtype_path, assay)
                for fname in os.listdir(complete_path):
                    fpath = complete_path + fname
                    if ' ' in fname:
                        new_fname = fname.replace(' ', '-')
                        fstem = fname.split('.')[0]
                        fext = fname.split('.')[1]
                        fstem = re.escape(fstem)
                        fname = "{}.{}".format(fstem,fext)
                        command = 'mv {}{} {}{}'.format(complete_path, fname, complete_path, new_fname)
                        print(command)
                        subprocess.call(command, shell=True)

if __name__=="__main__":
    fix_input_file_spaces()
