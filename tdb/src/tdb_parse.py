#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, re, datetime, gzip
from Bio import SeqIO
from Bio import Entrez
import requests
import pandas as pd
import numpy as np
from unidecode import unidecode

class tdb_parse(object):
    def __init__(self, **kwargs):
        pass

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
            sub_set_vals = HI_table[ref_serum][~np.isnan(HI_table[ref_serum])]
            sub_set_source = HI_table['source'][~np.isnan(HI_table[ref_serum])]
            sub_set_date = HI_table['collection date'][~np.isnan(HI_table[ref_serum])]
            sub_set_group = HI_table['genetic group'][~np.isnan(HI_table[ref_serum])]
            sub_set_passage = HI_table['passage history'][~np.isnan(HI_table[ref_serum])]
            sub_set_ref = HI_table['ref/test'][~np.isnan(HI_table[ref_serum])]
            for virus, val, src_id, date, passage, group, ref in zip(sub_set_vals.index, sub_set_vals, sub_set_source, sub_set_date, sub_set_passage, sub_set_group, sub_set_ref):
                flat_measurements.append({'virus': virus.upper(), 'serum': ref_serum[0].upper(), 'ferret_id': ref_serum[1].upper(), 'source': src_id.upper(), 'titer': val, 'date': date, 'passage': passage.upper(), 'group': group.upper(), 'ref': ref.upper(), 'date_modified': self.get_upload_date()})
        print("NIMR total:", len(flat_measurements), "measurements")
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
                        "NETH":"NETHERLANDS", "FIN":"FINLAND", "BRIS":"BRISBANE", "MARY":"MARYLAND",
                        "ST.P'BURG":"ST.PETERSBURG", 'CAL':'CALIFORNIA', 'CA':'CALIFORNIA', 'AUCK':'AUCKLAND', "C'CHURCH":'CHRISTCHURCH', "CCHURCH":'CHRISTCHURCH', 'CC': 'CHRISTCHURCH',
                        'CHCH':'CHRISTCHURCH', 'ASTR':'ASTRAKHAN', 'ASTRAK':'ASTRAKHAN', 'ST.P':"ST.PETERSBURG",'ST P':"ST.PETERSBURG",'STP':"ST.PETERSBURG", 'ST.PBURG': "ST.PETERSBURG", 'STPBURG':'ST.PETERSBURG',
                        'JHB':'JOHANNESBURG', 'FOR':'FORMOSA','MAL':'MALAYSIA', 'STHAUS':'SOUTHAUSTRALIA', 'BAY':'BAYERN',
                        'FL':'FLORIDA', 'MASS':'MASSACHUSETTS','NOVO':'NOVOSIBIRSK','WIS':'WISCONSIN','BANG':'BANGLADESH','EG':'EGYPT', 'SLOV':'SLOVENIA'}
        src_id = fname.split('/')[-1]
        with myopen(fname) as infile:
            csv_reader = csv.reader(infile)

            # parse sera
            row1 = csv_reader.next()
            row2 = csv_reader.next()
            row3 = csv_reader.next()
            if any('Passage history' in item for item in row3):
                row3 = csv_reader.next()
            ref_sera = [[HI_fix_name(e1+'/'+e2), e3.replace(' ', '')] for e1, e2, e3 in zip(row1, row2, row3)[4:]]
            for ri in xrange(len(ref_sera)):
                # replace abbreviations
                abbr = ref_sera[ri][0].split('/')[1].rstrip('01234566789')
                if abbr in name_abbrev:
                    ref_sera[ri][0] = HI_fix_name(ref_sera[ri][0].replace(abbr, name_abbrev[abbr]))
                else:
                    ref_sera[ri][0] = HI_fix_name(ref_sera[ri][0])
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
                    print(ref_sera[ri])

            fields = ['source','ref/test', 'genetic group', 'collection date', 'passage history']+map(tuple, ref_sera)
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
                    ref_strains.append(HI_fix_name(row[0].strip()))
                    ref_matrix.append([src_id,'ref']+map(strip, row[1:4])+map(titer_to_number, row[4:]))

            test_strains = []
            test_matrix = []
            for row in csv_reader: # load test viruses until it is no longer an A/ flu  name
                if not (row[0].startswith('A/') or row[0].startswith('B/')):
                    break
                else:
                    name = unidecode(row[0].strip().decode('utf-8'))
                    test_strains.append(HI_fix_name(name))
                    test_matrix.append([src_id,'test']+map(strip,row[1:4])+map(titer_to_number, row[4:]))
            HI_table  = pd.DataFrame(ref_matrix+test_matrix, index = ref_strains+test_strains, columns= fields)
            return HI_table

def titer_to_number(val):
    try:
        if '<' in val:
            return np.nan
        if len(val.split())>1:
            return float(val.split()[0])
        else:
            return float(val)
    except:
        #print "Bad HI measurement:", val
        return np.nan

def myopen(fname, mode='r'):
    if fname[-2:]=='gz':
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def HI_fix_name(name):
    if name.split() == ["NIB-85", "(A/Almaty/2958/2013)"]:
        print(name)
        tmp_name = fix_name("A/Almaty/2958/2013")
    elif name.split() == ["A/Texas/50/2012","(6&7)"]:
        tmp_name = fix_name("A/Texas/50/2012")
    elif name.split() == ["X-243", "(A/SOUTHAFRICA/3626/2013)"]:
        tmp_name = fix_name("A/SOUTHAFRICA/3626/2013")
    else:
        tmp_name = fix_name(name)
    return tmp_name.upper().lstrip('*')

def fix_name(name):
    tmp_name = name.replace(' ', '').replace('\'', '').replace('(', '').replace(')', '').replace('H3N2', '').replace('Human', '').replace('human', '').replace('//', '/').replace('.','')
    fields = tmp_name.split('/')
    if len(fields[-1])==2:
        try:
            y = int(fields[-1])
            if y>16:
                y=1900+y
            else:
                y=2000+y
            return '/'.join(fields[:-1])+'/'+str(y)
        except:
            return tmp_name
    else:
        return tmp_name