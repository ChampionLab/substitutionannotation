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
df['Multiple Sub'] = df['PTM'].str.count('substitution') >= 1
df['Number of Modifications'] = df['PTM'].str.count(';') - df['PTM'].str.count('Carbamidomethylation')

#Drop peptides without intensity in any sample
signalcols = df.columns[df.columns.str.contains('Area')]
nosignal = df.loc[:,signalcols].apply(lambda row: row.isna().all(),axis=1)
df = df[~nosignal]

"""Clarify Base and Mod peptide sequences"""
now = datetime.datetime.now()
print(str(now))
print('Clarify sequences')

#Get substitution origins, destinations, and positions
dfsubpositions = df['AScore'].str.extractall('(\d+):...->').reset_index(
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
#Also remove other, non substituted modification labels
df['Sub Peptide']=df['Peptide'].replace(regex=r'(.)\(sub [A-Z]\)',value=r'\1')
df['Sub Peptide'].replace(regex='\(.+?\)',value='',inplace = True)

#Extract base sequence from combined label in Peptide column, remove other modifications
df['Peptide'].replace(regex='(.\(sub (.)\))',value='\g<2>',inplace=True)
df['Peptide'].replace(regex='\(.+?\)',value='',inplace = True)


"""Output filtered PSM lists"""
now = datetime.datetime.now()
print(str(now))
print('Output Filtered PSMs')

#Rename/reshape df to unify downstream analysis

dict_col_names = {
    'RT':'Retention',
    'Length':'Peptide Length',
    'Source File':'Spectrum',
    'm/z':'Calibrated Observed M/Z',
    'Accession':'Protein',
    'Mass':'Calibrated Observed Mass',
    'Peptide':'Base Peptide',
    'Sub Peptide':'Peptide'
}
df.to_csv(os.path.join(sA.outputdir,'AllPSMSAndFilters.csv'))
#Get substitued peptide info
dfsubs=df[df['Is Sub']]

#Get substitution position in the protein

#Get quantification info
dfquant = dfsubs[['Protein','Sample','Peptide']]
dfsubs.to_csv(os.path.join(sA.outputdir,'SSP PSM.csv'))

dffilteredpsm[]


"""Output Total number of PSMs by mod, total sub PSMs by file"""

#Number of substitutions by file
bysample = df.groupby(by='Sample')
#Add back danger mods when implemented
dfsummary = bysample[['Is Sub','Intense Sub','Has Intensity']].sum()
dfsummary.loc[:,'Total'] = bysample['Is Sub'].size()
dfsummary.loc['Total',:]=dfsummary.sum()
dfsummary.to_csv(os.path.join(p.peaksout,'Mod PSM Summary by File.csv'))

#Number of PSMs by mod
"""The quick hack of getting most things"""
#spoiler alert it is not fast
modlist += ['Carbamidomethylation','Oxidation','Deamid']

now = datetime.datetime.now()
print(str(now))
print('Getting #PSMs by mod type and file')
"""Use loop for getting raw counts of PTMs"""

dfmodsum = pd.DataFrame(index=modlist,columns=filenames)
dfintensemodsum  = pd.DataFrame(index=modlist,columns=filenames)
for file in filenames:
    dfiter = df[df['Source File']==file]
    for mod in modlist:
        dfmodsum.loc[mod,file] = dfiter[dfiter['PTM'].str.contains(mod)]['Has Intensity'].count()
        dfintensemodsum.loc[mod,file] = dfiter[dfiter['PTM'].str.contains(mod)]['Has Intensity'].sum()

""""Or use groupby if you want all represented permutations of PTM combos. Currently untested"""
#byfilemod = df.groupby(by=['Source File','PTM'])['Has Intensity']
#dfmodsum = byfilemod.count()
#dfintensemodsum = byfilemod.sum()
#dfmodsum.reset_index(level=['Source File'])
#dfintensemodsum.reset_index(level=['Source File'])
dfmodsum.loc[:,'Total'] = dfmodsum.sum()
dfmodsum.loc['Total',:] = dfmodsum.sum(axis=1)
dfintensemodsum.loc[:,'Total'] = dfintensemodsum.sum()
dfintensemodsum.loc['Total',:] = dfintensemodsum.sum(axis=1)

dfmodsum.to_csv(os.path.join(p.peaksout,'Mod PSM Summary by Mod.csv'))
dfintensemodsum.to_csv(os.path.join(p.peaksout,'Intense Mod PSM Summary by Mod.csv'))

byfilemod = df.groupby(by=['Source File','PTM'])['Has Intensity']
dfmodsum = byfilemod.count()
dfintensemodsum = byfilemod.sum()





now = datetime.datetime.now()
print(str(now))
print('Done')

#Reset variables
#%reset -f