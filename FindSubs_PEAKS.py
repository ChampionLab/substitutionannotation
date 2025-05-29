#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:43:35 2021
Inputs - spyder.psms.csv
    
@author: taylorlundgren

To-do list:
    Add danger mod filtering...somehow?
    Expand ambiguous ->Leu substitution assignments

"""
"""Program to count the subs from PEAKS output"""

import pandas as pd
import os
import datetime
from substitutionannotation.resources import assets as sA
import numpy as np

#Create list of every possible amino acid substitution
aminoacids = ['Ala', 'Arg','Asn','Asp','Cys','Glu','Gln','Gly',
              'His','Lys','Met','Phe','Pro','Ser',
              'Thr','Trp','Tyr','Val','Leu','Ile']
dict_renamePos = { 0: 'Sub 1 Position',
               1 : 'Sub 2 Position',
               2 : 'Sub 3 Position',
               3 : 'Sub 4 Position',
               4 : 'Sub 5 Position'
               }
dict_renameOri = { 0: 'Sub 1 Origin',
               1 : 'Sub 2 Origin',
               2 : 'Sub 3 Origin',
               3 : 'Sub 4 Origin',
               4 : 'Sub 5 Origin',
               }
dict_renameDest = { 0: 'Sub 1 Destination',
               1 : 'Sub 2 Destination',
               2 : 'Sub 3 Destination',
               3 : 'Sub 4 Destination',
                4: 'Sub 5 Destination'
               }

modlist=[]
for i in aminoacids:
    for j in aminoacids:
        if i!=j:
            modlist.append(i +'->'+ j)


"""-----------------------------Import-------------------------------------"""
now = datetime.datetime.now()
print(str(now))
print('Start import')

df = pd.read_csv(os.path.join(sA.inputdir,'spider.peptides.csv'))
samples = pd.read_csv(os.path.join(sA.inputdir,'spider.filteredStatistics.csv'))['Sample'].iloc[:-1]

#df['Full Spectrum'] = [str(x) + '-' + y for x, y in zip(df['Scan'], df['Source File'])]

df.loc[:,'PTM'].fillna('NONE',inplace=True)
#df['Sample'] = df['Source File'].replace(regex='(?=_Slot)(.+)',value='')
#filenames = df['Source File'].replace(regex='(?=_Slot)(.+)',value='').drop_duplicates()
#filenames.to_csv(os.path.join(p.peaksout,'filenames.csv'),index=False)

now = datetime.datetime.now()
print(str(now))
print('Start filters')

df['Is Sub'] = df['PTM'].str.count('substitution') == 1
df['Multiple Sub'] = df['PTM'].str.count('substitution') > 1
df['Number of Modifications'] = df['PTM'].str.count(';') - df['PTM'].str.count('Carbamidomethylation') + 1

#Drop peptides without intensity in any sample
signalcols = df.columns[df.columns.str.contains('Area')]
nosignal = df.loc[:,signalcols].apply(lambda row: row.isna().all(),axis=1)
df = df[~nosignal]

"""Clarify Base and Mod peptide sequences"""
now = datetime.datetime.now()
print(str(now))
print('Clarify sequences')

#Get substitution origins, destinations, and positions
dfsubpositions = df['AScore'].str.extractall(r'(\d+):...->').reset_index(
    level='match').pivot(columns ='match', values=0).rename(
        columns=dict_renamePos)
df = pd.concat([df,dfsubpositions],axis=1)

dfsuborigins = df['AScore'].str.extractall('(...)->').replace(sA.dict_321
        ).reset_index(level='match').pivot(columns ='match', values=0).rename(
        columns=dict_renameOri)
df = pd.concat([df,dfsuborigins],axis=1)

dfsubdestinations = df['AScore'].str.extractall('->(...)').replace(sA.dict_321
        ).reset_index(level='match').pivot(columns ='match', values=0).rename(
        columns=dict_renameDest)
df = pd.concat([df,dfsubdestinations],axis=1)

#Extract substituted sequence from combined label in Peptide column
df['Sub Peptide']=df['Peptide'].replace(regex=r'(.)\(sub [A-Z]\)',value=r'\1')

#Also remove other, non substituted modification labels
df['Sub Peptide'].replace(regex=r'\(.+?\)',value='',inplace = True)

#Extract base sequence from combined label in Peptide column, remove other modifications
df['Base Peptide'] = df['Peptide'].replace(regex=r'(.\(sub (.)\))',value=r'\g<2>')
df['Base Peptide'].replace(regex=r'\(.+?\)',value='',inplace = True)


"""Output filtered PSM lists"""
now = datetime.datetime.now()
print(str(now))
print('Output Filtered PSMs')

#Rename/reshape df to unify downstream analyses

dict_col_names = {
    'RT':'Retention',
    'Length':'Peptide Length',
    'Source File':'Spectrum',
    'm/z':'Calibrated Observed M/Z',
    'Accession':'Protein',
    'Mass':'Calibrated Observed Mass',
    'Sub Peptide':'Modified Peptide'
}
df.rename(columns=dict_col_names,inplace=True)
df.to_csv(os.path.join(sA.outputdir,'AllPSMSAndFilters.csv'))
#Get substitued peptide info
dfsubs=df[df['Is Sub']]

#Get substitution position in the protein

#Format quantification info
id_cols=['Protein','Peptide','Modified Peptide','Base Peptide']
area_cols = [f'Area {sample}' for sample in samples]
dfquant = dfsubs.melt(id_vars=id_cols,value_vars=area_cols,var_name='Sample',value_name='Intensity')
dfquant['Sample'].replace(regex='Area ',value='',inplace=True)
#Get base peptide info
dfbase = df[df['Peptide'].isin(dfquant['Base Peptide'].unique())].melt(id_vars=['Peptide'],value_vars=area_cols,var_name='Sample',value_name='Base Intensity')
dfbase.rename(columns={'Peptide':'Base Peptide'},inplace=True)
dfbase['Sample'].replace(regex='Area ',value='',inplace=True)

dfquant = dfquant.merge(dfbase,on=['Base Peptide','Sample'], how='left')
#Calculate Values
dfquant['Ratio'] = dfquant['Intensity']/dfquant['Base Intensity']
dfquant['logRatio'] = np.log10(dfquant['Ratio'].replace(0,np.nan))
dfsubs.to_csv(os.path.join(sA.outputdir,'SSP PSM.csv'))
dfquant.to_csv(os.path.join(sA.outputdir,'SSP Quant.csv'))



now = datetime.datetime.now()
print(str(now))
print('Done')

#Reset variables
#%reset -f