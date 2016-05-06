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

class tdb_parse(object):
    def __init__(self, **kwargs):
        self.table_column_names = ['viruses', 'other', 'collection', 'passage', '']
        self.titer_values = ['10.0', '20.0','40.0', '80.0', '160.0', '320.0', '640.0', '1280.0', '2560.0', '5120.0', 'nan']

    def parse(self):
        '''
        Parse HI data files to create list of measurement documents
        '''
        HI_titers = self.read_tables()
        self.measurements = self.table_to_flat(HI_titers)

    def get_upload_date(self):
        '''
        return current date
        '''
        return str(datetime.datetime.strftime(datetime.datetime.now(),'%Y-%m-%d'))

    def read_tables(self):
        '''
        Read all csv tables in path, create data frame with reference viruses as columns
        '''
        import glob
        flist = glob.glob(self.path + '/NIMR*csv')
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
                    flat_measurements.append({'virus': virus.upper(), 'serum': ref_serum[0].upper(), 'ferret_id': ref_serum[1].upper(), 'source': src_id.upper(), 'titer': val, 'date': date, 'passage': passage.upper(), 'ref': ref.upper(), 'date_modified': self.get_upload_date()})
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
        name_abbrev = {'HK':"HONGKONG", 'SWITZ':"SWITZERLAND", 'VIC':"VICTORIA", 'STOCK':"STOCKHOLM", 'DR':'DOMINICANREPUBLIC', 'NJ': 'NEWJERSEY',
                        'STHAFR':"SOUTHAFRICA", 'SAFRICA':"SOUTHAFRICA", "ENG":"ENGLAND", "NIB-85":"A/ALMATY/2958/2013", 'X-243':'A/SOUTHAFRICA/3626/2013', 'NOR':'NORWAY',
                        'NTHCAROL':"NORTHCAROLINA",'ALA':"ALABAMA", 'NY':"NEWYORK", "GLAS":"GLASGOW", "AL":"ALABAMA",
                        "NETH":"NETHERLANDS", "FIN":"FINLAND", "BRIS":"BRISBANE", "MARY":"MARYLAND", "BRIS1,": "BRISBANE", "MAD": "MADAGASCAR", "SAUSTRALIA": "SOUTHAUSTRALIA", "SHAN" : "SHANDONG",
                        "ST.P'BURG":"ST.PETERSBURG", 'CAL':'CALIFORNIA', 'CA':'CALIFORNIA', 'AUCK':'AUCKLAND', "C'CHURCH":'CHRISTCHURCH', "CCHURCH":'CHRISTCHURCH', 'CC': 'CHRISTCHURCH',
                        'CHCH':'CHRISTCHURCH', 'ASTR':'ASTRAKHAN', 'ASTRAK':'ASTRAKHAN', 'ST.P':"ST.PETERSBURG",'ST P':"ST.PETERSBURG",'STP':"ST.PETERSBURG", 'ST.PBURG': "ST.PETERSBURG", 'STPBURG':'ST.PETERSBURG',
                        'JHB':'JOHANNESBURG', 'FOR':'FORMOSA','MAL':'MALAYSIA', 'STHAUS':'SOUTHAUSTRALIA', 'BAY':'BAYERN', "BAN": "BANGLADESH", "VALL": "VALLADOLID", "NIEDER":"NIEDERSACHSEN",
                        'FL':'FLORIDA', 'MASS':'MASSACHUSETTS','NOVO':'NOVOSIBIRSK','WIS':'WISCONSIN','BANG':'BANGLADESH','EG':'EGYPT', 'SLOV':'SLOVENIA',
                        "H-X": "Heilongjiang-Xiangfang", "VORO": "Voronezh"}
        src_id = fname.split('/')[-1]
        with self.myopen(fname) as infile:
            csv_reader = csv.reader(infile)
            # parse sera
            row1 = csv_reader.next()
            row2 = csv_reader.next()
            row3 = csv_reader.next()
            fields = self.determine_table_meta_fields(row1)
            ref_sera_start = len(fields)
            if not all(field in fields for field in ['collection', 'passage']):
                fields = self.determine_table_meta_fields(row2)
                fields[0] = 'viruses'
            ref_sera_start = len(fields)
            try:
                fields.remove("viruses")
            except:
                print("couldn't remove viruses field", fields, src_id)
            # starting 2016, included passage history row
            if re.match(r'(\D+)(\d\d\d\d)', src_id):  #NIMR_Feb2012_10.csv, NIMR-report-Feb2011_04.csv
                num = re.match(r'(\D+)(\d\d\d\d)', src_id).group(2)
                try:
                    num = int(num)
                    if num >= 2016:
                        row3 = csv_reader.next()
                except:
                    pass
                    print("couldn't parse file name")
            row3 = [re.match(r'^([^\*]*)', id).group(0).upper() for id in row3]  # get everything before the '*'?
            ref_sera = [[self.HI_fix_name(e1+'/'+e2), e3.replace(' ', '')] for e1, e2, e3 in zip(row1, row2, row3)[ref_sera_start:]]
            for ri in xrange(len(ref_sera)):
                # replace abbreviations
                abbr = ref_sera[ri][0].split('/')[1].rstrip('01234566789')
                if abbr in name_abbrev:
                    ref_sera[ri][0] = self.HI_fix_name(ref_sera[ri][0].replace(abbr, name_abbrev[abbr]))
                else:
                    ref_sera[ri][0] = self.HI_fix_name(ref_sera[ri][0])
                # strip numbers
                tmp = ref_sera[ri][0].split('/')
                ref_sera[ri][0] = '/'.join([tmp[0], tmp[1].rstrip('0123456789')]+tmp[2:])
                try:
                    y = int(ref_sera[ri][0].split('/')[-1])
                    if y<100:
                        if y<20:
                            ref_sera[ri][0] = '/'.join(ref_sera[ri][0].split('/')[:-1])+'/'+str(2000+y)
                        else:
                            ref_sera[ri][0] = '/'.join(ref_sera[ri][0].split('/')[:-1])+'/'+str(1900+y)
                except:
                    print("Couldn't parse this reference name", ref_sera[ri], src_id)
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
                    ref_strains.append(self.HI_fix_name(row[0].strip()))
                    ref_matrix.append([src_id,'ref']+map(strip, row[1:ref_sera_start])+map(self.titer_to_number, row[ref_sera_start:]))

            test_strains = []
            test_matrix = []
            for row in csv_reader:
                try:
                    name = unidecode(row[0].strip().decode('utf-8'))
                    test_strains.append(self.HI_fix_name(name))
                    test_matrix.append([src_id,'test']+map(strip,row[1:ref_sera_start])+map(self.titer_to_number, row[ref_sera_start:]))
                    # look for titer values that are not normal
                    titers = map(self.titer_to_number, row[ref_sera_start:])
                    for t in titers:
                        t = str(t)
                        if t not in self.titer_values:
                            print("Weird titer value", t, src_id)
                except:
                    print("Couldn't parse name from file", row[0].strip(), src_id)
            HI_table  = pd.DataFrame(ref_matrix+test_matrix, index = ref_strains+test_strains, columns= fields)
            # get rid of columns 'other' and ''
            if 'other' in HI_table:
                HI_table = HI_table.drop('other', 1)
            while '' in HI_table:
                HI_table = HI_table.drop('', 1)
            return HI_table
    def determine_table_meta_fields(self, row1):
        fields = []
        for col in row1:
            if col.strip().lower() not in self.table_column_names:
                break
            else:
                fields.append(col.strip().lower())
        return fields

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

    def myopen(self, fname, mode='rU'):
        if fname[-2:]=='gz':
            return gzip.open(fname, mode)
        else:
            return open(fname, mode)

    def HI_fix_name(self, name):
        if name.split() == ["NIB-85/"]:
            tmp_name = self.fix_name("A/Almaty/2958/2013")
        elif name.split() == ["A/Texas/50/2012","(6&7)"]:
            tmp_name = self.fix_name("A/Texas/50/2012")
        elif name.split() == ["X-243/"]:
            tmp_name = self.fix_name("X-243A/SOUTHAFRICA/3626/2013")
        else:
            tmp_name = self.fix_name(name)
        return tmp_name.upper().lstrip('*')

    def fix_name(self, name):
        tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.', '').replace(',', '')
        fields = tmp_name.split('/')
        lookup_month = {'Jan': '1', 'Feb': '2', 'Mar': '3', 'Apr': '4', 'May': '5', 'Jun': '6',
                            'Jul': '7', 'Aug': '8', 'Sep': '9', 'Oct': '10', 'Nov': '11', 'Dec': '12'}
        try:
            if len(fields) ==4:
                fields[2] = fields[2].lstrip('0')
                tmp_name = '/'.join(fields)
        except:
            pass
        if re.match(r'(\d+)-(\D+)',fields[-1]):  # A/CALIFORNIA/9-APR
            try:
                month = lookup_month[re.match(r'(\d+)-(\D+)',fields[-1]).group(2).title()]
                year = re.match(r'(\d+)-(\D+)',fields[-1]).group(1)
                if len(year) == 1:
                    y = '200' + year
                elif len(year) == 2:
                    y = int(year)
                    if y>16:
                        y=1900+y
                    else:
                        y=2000+y
                return '/'.join(fields[:-1]) + '/' + month + '/' + str(y)
            except:
                return tmp_name
        elif re.match(r'(\D+)-(\d+)',fields[-1]):  # B/SHANDONG/JUL-97
            try:
                month = lookup_month[re.match(r'(\D+)-(\d+)',fields[-1]).group(1).title()]
                year = re.match(r'(\D+)-(\d+)',fields[-1]).group(2)
                if len(year) == 1:
                    y = '200' + year
                elif len(year) == 2:
                    y = int(year)
                    if y>16:
                        y=1900+y
                    else:
                        y=2000+y
                return '/'.join(fields[:-1]) + '/' + month + '/' + str(y)
            except:
                return tmp_name
        elif len(fields[-1])==2:
            try:
                y = int(fields[-1])
                if y>16:
                    y=1900+y
                else:
                    y=2000+y
                return '/'.join(fields[:-1])+'/'+str(y)
            except:
                return tmp_name
        elif re.match(r'[\d,\*\D\-]+([AB]$)', fields[0]): # 1,3B/FLORIDA/4/2006  #NYMCX-263BA/HK/4801/2014
            try:
                result = re.match(r'[\d,\*\D\-]+([AB]$)', fields[0]).group(1)
                return result + '/' + '/'.join(fields[1:])
            except:
                return tmp_name
        elif re.match(r"[\s\w\-]+\(([/\w]+)\)", tmp_name): # IVR-159 (A/Victoria/502/2010)
            try:
                return re.match(r"[\s\w\-]+\(([/\w]+)\)", tmp_name).group(1)
            except:
                return tmp_name
        else:
            return tmp_name