#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, re, datetime, gzip
from Bio import SeqIO
from Bio import Entrez
import requests
import pandas as pd
import numpy as np
from unidecode import unidecode
import math

class parse(object):
    def __init__(self, **kwargs):
        self.table_column_names = ['viruses', 'other', 'collection', 'passage', '']
        self.titer_values = ['10.0', '20.0', '40.0', '80.0', '160.0', '320.0', '640.0', '1280.0', '2560.0', '5120.0', 'nan']

    def parse(self, ftype='flat', **kwargs):
        '''
        Parse HI data files to create list of measurement documents
        '''
        flat_measurements = list()
        if ftype == "flat":
            flat_measurements = self.read_flat(**kwargs)
        elif ftype == "tables":
            HI_titers = self.read_tables(**kwargs)  
            flat_measurements = self.table_to_flat(HI_titers)
        return flat_measurements
        
    def read_flat(self, path, fstem, **kwargs):
        '''
        Read flat titer table, assumes file ends with .tsv
        Example line:
        A/SYDNEY/5/1997	A/NETHERLANDS/22/2003	NL/22/03	Smith2004	320.0
        '''
        import csv
        file = path + fstem + ".tsv"
        flat_measurements = list()
        try:
            os.path.isfile(file)
        except IOError:
            raise Exception(file, "not found")
        else:
            with open(file) as infile:
                table_reader = csv.reader(infile, delimiter="\t")
                header = {
                    0: 'virus_strain',
                    1: 'serum_strain',
                    2: 'ferret_id',
                    3: 'source',
                    4: 'titer'
                }
                for row in table_reader:
                    m = {key: row[ii] if ii < len(row) else "" for ii, key in header.items()}
                    m['passage'] = 'unknown'  # TODO FIX TOTAL HACK
                    flat_measurements.append(m)
        return flat_measurements        

    def read_tables(self, path, **kwargs):
        '''
        Read all csv tables in path, create data frame with reference viruses as columns
        '''
        import glob
        flist = glob.glob(path + '/NIMR*csv')
        HI_matrices = pd.DataFrame()
        for fname in flist:
            tmp = self.parse_HI_matrix(fname)
            HI_matrices = HI_matrices.append(tmp)
        return HI_matrices

    def table_to_flat(self, HI_table):
        flat_measurements = list()
        for ref_serum in HI_table.columns[5:]:
            try:
                sub_set_vals = HI_table[ref_serum][~np.isnan(HI_table[ref_serum])]
                sub_set_source = HI_table['source'][~np.isnan(HI_table[ref_serum])]
                sub_set_date = HI_table['collection'][~np.isnan(HI_table[ref_serum])]
                sub_set_passage = HI_table['passage'][~np.isnan(HI_table[ref_serum])]
                sub_set_ref = HI_table['ref/test'][~np.isnan(HI_table[ref_serum])]
                for virus, val, src_id, date, passage, ref in zip(sub_set_vals.index, sub_set_vals, sub_set_source, sub_set_date, sub_set_passage, sub_set_ref):
                    flat_measurements.append({'virus_strain': virus, 'serum_strain': ref_serum[0], 'ferret_id': ref_serum[1], 'source': src_id, 'titer': val, 'date': date, 'passage': passage, 'ref': ref})
            except:
                print("Couldn't parse this serum's measurements", ref_serum)
                print("Check fields at top left of file")
        return flat_measurements

    def parse_HI_matrix(self, fname):
        '''
        Parse HI file to create dataframe with reference viruses as columns
        :param fname:
        :return:
        '''
        from string import strip
        import csv
        src_id = fname.split('/')[-1]
        with self.myopen(fname) as infile:
            csv_reader = csv.reader(infile)
            # parse sera
            row1 = csv_reader.next()
            row2 = csv_reader.next()
            row3 = csv_reader.next()
            # starting 2016, included passage history row, need to get fourth row
            if self.determine_source_year(src_id) >= 2016:
                row3 = csv_reader.next()
            fields = self.determine_columns(row1)
            ref_sera_start = len(fields)
            if not all(field in fields for field in ['collection', 'passage']):
                fields = self.determine_columns(row2)
                fields[0] = 'viruses'
            ref_sera_start = len(fields)
            try:
                fields.remove("viruses")
            except:
                print("couldn't remove viruses field", fields, src_id)
            row3 = [re.match(r'^([^\*]*)', id).group(0).upper() for id in row3]  # get everything before the '*'?
            ref_sera = [[(e1+'/'+e2), e3.replace(' ', '')] for e1, e2, e3 in zip(row1, row2, row3)[ref_sera_start:]]
            fields = ['source','ref/test'] + fields + map(tuple, ref_sera)
            #print fields
            for row in csv_reader: # advance until the reference virus
                if row[0].startswith('REFERENCE'):
                    break
            ref_strains = []
            ref_matrix = []
            for row in csv_reader:
                if row[0].startswith('TEST'):
                    break
                else: # load matrices until the test virus section starts
                    ref_strains.append(row[0].strip())
                    ref_matrix.append([src_id,'ref']+map(strip, row[1:ref_sera_start])+map(self.titer_to_number, row[ref_sera_start:]))

            test_strains = []
            test_matrix = []
            for row in csv_reader:
                try:
                    name = unidecode(row[0].strip().decode('utf-8'))
                    test_strains.append(name)
                    test_matrix.append([src_id,'test']+map(strip,row[1:ref_sera_start])+map(self.titer_to_number, row[ref_sera_start:]))
                    self.check_titer_values(map(self.titer_to_number, row[ref_sera_start:]), src_id)
                except:
                    print("Couldn't parse name from file", row[0].strip(), src_id)

            HI_table  = pd.DataFrame(ref_matrix+test_matrix, index = ref_strains+test_strains, columns= fields)
            # get rid of columns 'other' and ''
            if 'other' in HI_table:
                HI_table = HI_table.drop('other', 1)
            while '' in HI_table:
                HI_table = HI_table.drop('', 1)
            return HI_table

    def determine_columns(self, row1):
        fields = []
        for col in row1:
            if col.strip().lower() not in self.table_column_names:
                break
            else:
                fields.append(col.strip().lower())
        return fields

    def determine_source_year(self, src_id):
        '''
        # starting 2016, included passage history row
        '''
        year = 0
        if re.match(r'\D+(\d\d\d\d)', src_id):  #NIMR_Feb2012_10.csv, NIMR-report-Feb2011_04.csv
            year = re.match(r'\D+(\d\d\d\d)', src_id).group(1)
            try:
                year = int(year)
            except:
                print("couldn't source file name to get year")
        return year

    def titer_to_number(self, val):
        try:
            if '<' in val:
                return np.nan
            elif '>' in val:
                return float(val)
            elif re.match(r'0\s+([0-9]+)', val): # 0 160
                val = re.match(r'0\s+([0-9]+)', val).group(1)
                return float(val)
            elif val + '.0' not in self.titer_values:
                temp = str(float(val) * 10.0)
                if temp in self.titer_values:
                    return float(temp)
            return float(val)
        except:
            #print "Bad HI measurement:", val
            return np.nan

    def check_titer_values(self, titers, src_id):
        '''
        look for titer values that are not normal
        '''
        for t in titers:
            t = str(t)
            if t not in self.titer_values:
                print("Weird titer value", t, src_id)

    def myopen(self, fname, mode='rU'):
        if fname[-2:]=='gz':
            return gzip.open(fname, mode)
        else:
            return open(fname, mode)
