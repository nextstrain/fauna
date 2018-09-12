import os
import sys
import xlrd
import re

output_file = "../source-data/vidrl_serum_mapping.tsv"

vidrl_excel_directory = "../../../fludata/VIDRL-Melbourne-WHO-CC/raw-data/"

sera_mapping = {}

for ab in os.listdir(vidrl_excel_directory):
        ab_path = '{}{}/'.format(vidrl_excel_directory, ab)
        for subtype in os.listdir(ab_path):
            subtype_path = '{}{}/'.format(ab_path, subtype)
            for assay in os.listdir(subtype_path):
                complete_path = '{}{}/'.format(subtype_path, assay)
                for fname in os.listdir(complete_path):
                    print("Parsing {}".format(fname))
                    if ' ' in fname:
                        print('\tbad file name, passing')
                        pass
                    fpath = complete_path + fname
                    workbook = xlrd.open_workbook(fpath)
                    for sheet in workbook.sheets():
                        try:
                            temp = sheet.row_values(10)
                            col = sheet.col_values(2)
                            for i in range(4,15):
                                t = temp[i].strip('\n')
                                if (type(t) != float):
                                    sera_mapping[t] = col[i+8]
                            print('\tSuccess!')
                        except:
                            print('\t\tfailed')

print(sera_mapping)

with open(output_file, 'w') as f:
    for key in sera_mapping.keys():
        f.write('{}\t{}\n'.format(key, sera_mapping[key]))
