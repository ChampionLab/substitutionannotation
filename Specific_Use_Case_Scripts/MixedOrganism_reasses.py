#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2024/02/26

Find and report on which singly substituted peptides (SSP) arising from the mixed
organism experiment were found and compare against the DB search
Compacted script to repeat some of the analyses done on the E. coli / S. typhimurium experiment on more recent MSFragger searches
@author: taylorlundgren
"""
#%%
from substitutionannotation.resources import assets as p
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  
import datetime
import os
import seaborn as sns
import re
import plotly.graph_objects as go
import ast
inputdir = input('Input directory:')
outputdir = os.path.join(inputdir,'Discovery_v_DB UpdatedDB')
if not os.path.isdir(outputdir):
    os.mkdir(outputdir)

dict_dilutions = {'2021-06-18-SALTY':'100',
                  '2021-06-18-ECLandSALTY1':'150',
                  '2021-06-18-ECLandSALTY2':'75',
                  '2021-06-18-ECLandSALTY3':'38',
                  '2021-06-18-ECLandSALTY3-5':'19',
                  '2021-06-18-ECLandSALTY3_5':'19',
                  '2021-06-18-ECLandSALTY4':'10',
                  '2021-06-18-ECLandSALTY5':'5',
                  '2021-06-18-ECLandSALTY6':'2.5',
                  '2021-06-18-ECLandSALTY7':'1.25',
                  '2021-06-18-ECLandSALTY8':'0.68',
                  '2021-06-18-ECLandSALTY9':'0.32',
                  '2021-06-18-ECL':'0.5'
                      }

dict_db_to_sample = {'2021-06-18-SALTY':100,
                  '2021-06-18-ECLandSALTY1':66,
                  '2021-06-18-ECLandSALTY2':50,
                  '2021-06-18-ECLandSALTY3':33,
                  '2021-06-18-ECLandSALTY3-5':20,
                  '2021-06-18-ECLandSALTY3_5':20,
                  '2021-06-18-ECLandSALTY4':11,
                  '2021-06-18-ECLandSALTY5':6,
                  '2021-06-18-ECLandSALTY6':3,
                  '2021-06-18-ECLandSALTY7':2,
                  '2021-06-18-ECLandSALTY8':1,
                  '2021-06-18-ECLandSALTY9':0.5,
                  '2021-06-18-ECL':0
                      }

dict_update_samples = {'SALTY':100, '66':66,'50':50,'33':33,'20':20,'11':11,'6':6,'3':3,'2':2,'1':1,
                          '0':0.5,
                  'ECL':0, 'ECOLI':0
                      }

order = ['100','66','50','33','20','11','6','3','2','1','0.5','0']
samples = [100,66,50,33,20,11,6,3,2,1,0.5,0]

#%%


print('Starting imports...')
now = datetime.datetime.now()
print(str(now))
""""Imports """
#substituted PSMs from search data
dfsubs = pd.read_csv(os.path.join(inputdir,'SSP PSM.csv'))
dfsubs['Sample'] = dfsubs['Sample'].astype(str)
dfsubs['Sample'].replace(dict_update_samples,inplace=True)
dfquant = pd.read_csv(os.path.join(inputdir,'SSP Quant.csv'))
dfquant['Sample'] = dfquant['Sample'].astype(str)
dfquant['Sample'].replace(dict_update_samples,inplace=True)
#All psms from discovery search
allpsm = pd.read_csv(os.path.join(inputdir,'AllPSMsandFilters.csv'),low_memory=False)
allpsm['Sample'] = allpsm['Sample'].astype(str)
allpsm['Sample'].replace(dict_update_samples,inplace=True)

#Add substitution position to allpsm
def find_sub_pos(pep):
    if type(pep) == float:
        return None
    if pep:
        res = re.search('.\[',pep)
        if res:
            return res.start()

allpsm['Substitution Peptide Position'] =  allpsm['Modified Peptide'].apply(
    find_sub_pos)

samplenames = dfsubs['Sample'].drop_duplicates()

#Orient target list based on provided FASTA
dfsummary = pd.DataFrame(index=samplenames)
if p.dbused =='SALTY':
    dftargets = p.targetlist[~p.targetlist['Is Exact']]
    dbfiltercol = 'ECOLI Pep Rep SSP'
elif p.dbused =='ECOLI':
    dftargets = p.targetlistrev[~p.targetlistrev['Is Exact']]
    dbfiltercol = 'SALTY Pep Rep SSP'
else:
    raise Exception('Please check the p.dbused variable to be SALTY or ECOLI')

#Pandas doesn't import lists very well...
dftargets['SSP Positions'] = dftargets['SSP Positions'].apply(ast.literal_eval)
dftargets['Substitution Types'] = dftargets['Substitution Types'].apply(ast.literal_eval)


#%%
#Determine quantification rate of substitutions
bySampleMod = dfsubs.groupby(by=['Sample','Modified Peptide'])
byMod = dfsubs.groupby(by=['Modified Peptide'])
bySample = dfsubs.groupby(by=['Sample'])
quantrateSample = bySample['Intensity'].apply(lambda x: (x>0).sum())/bySample['Intensity'].apply(len)
quantrate = (dfsubs['Intensity']>0).sum()/ len(dfsubs)
print(quantrateSample)
print(f'Global quantification rate: {quantrate}')
plt.bar(quantrateSample.index.astype(str),quantrateSample)
#plt.hlines(quantrate,0,11,'k','--','Global average')
plt.ylim(0,1)
plt.ylabel('% PSM Quantified')
plt.xlabel('% SALTY')
plt.savefig(os.path.join(outputdir,'Quantification Efficiency.png'))

#%%
#import concatenated DB closed search
print('Getting DB search results...')
now = datetime.datetime.now()
print(str(now))

columns = ['Peptide','Spectrum','Intensity','Hyperscore', 'Number of Missed Cleavages',
           'Probability','Calibrated Observed M/Z'
           ]
#Get DB search file/sample names
# for old DB search filenames = pd.read_csv(r"\\mchampion-nas.esc.nd.edu\DATA\TJL\Python Backups\Analysis Backups\MixedExperiment\Final\DefaultDB\filelist_ionquant.txt",sep='\t')
filenames = pd.read_csv(r"\\mchampion-nas.esc.nd.edu\DATA\TJL\TIMS\2021-06-18\MSFragger_Searches\v22\CombinedDB\filelist_ionquant.txt",sep='\t')
filenames = filenames[filenames['flag'] == '--psm']
filenames = filenames['value'].apply(os.path.dirname)
#Concatenate DB search PSM files
dfdb = pd.DataFrame()
columns += ['Retention','Ion Mobility','Protein']
dbsearchdir = r'\\mchampion-nas.esc.nd.edu\DATA\TJL\TIMS\2021-06-18\MSFragger_Searches\v22\CombinedDB'
for file in filenames:
    try: dfiter = pd.read_csv(os.path.join(dbsearchdir,file,'psm.tsv'), sep='\t', usecols=columns)
    except: print(file+" did not import properly")
    else: dfdb = pd.concat([dfdb,dfiter])
dfdb.reset_index(drop=True,inplace=True)

#Sometimes fragger hangs up on filenames with a .
#For runs where this needs to be renamed, I need to rename the spectrum data as well.
dfdb['Spectrum'] = dfdb['Spectrum'].str.replace(r'SALTY3.5_','SALTY3_5_')
dfdb['Spectrum'] = dfdb['Spectrum'].str.replace(r'SALTY3-5_','SALTY3_5_')

#Identify target SAAVs
dfdb.loc[:,'Is Target'] = dfdb['Peptide'].isin(dftargets['Base Sequence'])
#dfdb = dfdb[(dfdb['SALTY Pep Rep SSP']|dfdb['ECOLI Pep Rep SSP'])]
dfdb['Sample'] = dfdb['Spectrum'].replace(regex='(?=_Slot)(.+)',value='')
dfdb['Sample'] = dfdb['Sample'].replace(dict_db_to_sample)
dfdb['Sample'] = dfdb['Sample'].replace(dict_update_samples)    #This second dictionary renaming fixes updated DB search
alldb = dfdb.copy()
dfdb = alldb[alldb['Is Target']]
dfdb.columns = 'DB ' + dfdb.columns

#%%
#Merge DB and sub info

merged = dfdb.merge(allpsm,how='left',left_on='DB Spectrum', right_on='Spectrum')
shared = merged.copy()
targetsbybase = dftargets.groupby(by='Base Sequence')
#Categorize mismatches
    #Note that the inputs have duplicated values for I/L ambiguities...
noid = merged[merged['Is Sub'].isna()]          #Spectra without any ID in discovery
merged = merged[~merged['Is Sub'].isna()]

def good_base(row):
    return row['Peptide'] in targetsbybase.get_group(row['DB Peptide'])['ECOLI SSP Sequence'].iloc[0]
good_base_bool = merged.apply(good_base,axis=1)
bad_base = merged[~good_base_bool]                #Spectra with incorrect unmodified sequence
merged = merged[good_base_bool]

no_mod = merged[merged['PTMs'] == 'NONE']       #Spectra which were left unmodified
merged = merged[~(merged['PTMs'] == 'NONE')]

other_mod = merged[~merged['PTMs'].str.contains('->')]          #Spectra ONLY assigned a non-substitution or unknown modification
merged = merged[merged['PTMs'].str.contains('->')]

ambiguous = merged[merged['PTMs'].str.contains('or')]           #Spectra whose substitution was ambiguous with anoter modification
merged = merged[~merged['PTMs'].str.contains('or')]

def good_loc(row):
    possible_positions = targetsbybase.get_group(row['DB Peptide'])['SSP Positions'].iloc[0]
    return row['Substitution Peptide Position'] in possible_positions
good_loc_bool = merged.apply(good_loc,axis=1)
bad_loc = merged[~good_loc_bool]                                #Spectra with incorrectly localized substitution
merged = merged[good_loc_bool]


wrong_sub = merged[~(merged['DB Peptide'] == merged['Substituted Sequence'])]          #Spectra correctly assigned a substitution, but determined an incorrect sequence
merged = merged[merged['DB Peptide'] == merged['Substituted Sequence']]

no_int = merged[~(merged['Intensity'] > 0)]                     #Substitutions that got the right sequence, but no intensity
correct = merged[merged['Intensity'] > 0]  

#Handle duplication of rows for I/L ambiguities through hierarchy of categorization
    #Remove duplicates within each subsection
shared = shared.drop_duplicates(subset='DB Spectrum')
correct = correct.drop_duplicates(subset='DB Spectrum')
no_int = no_int.drop_duplicates(subset='DB Spectrum')
wrong_sub = wrong_sub.drop_duplicates(subset='DB Spectrum')
bad_loc = bad_loc.drop_duplicates(subset='DB Spectrum')
ambiguous = ambiguous.drop_duplicates(subset='DB Spectrum')
other_mod = other_mod.drop_duplicates(subset='DB Spectrum')
no_mod = no_mod.drop_duplicates(subset='DB Spectrum')
bad_base = bad_base.drop_duplicates(subset='DB Spectrum')
noid = noid.drop_duplicates(subset='DB Spectrum')
    #Remove duplicates found in correct spectra
no_int = no_int[~no_int['DB Spectrum'].isin(correct['DB Spectrum'])]
wrong_sub = wrong_sub[~wrong_sub['DB Spectrum'].isin(correct['DB Spectrum'])]
bad_loc = bad_loc[~bad_loc['DB Spectrum'].isin(correct['DB Spectrum'])]
ambiguous = ambiguous[~ambiguous['DB Spectrum'].isin(correct['DB Spectrum'])]
other_mod = other_mod[~other_mod['DB Spectrum'].isin(correct['DB Spectrum'])]
no_mod = no_mod[~no_mod['DB Spectrum'].isin(correct['DB Spectrum'])]
bad_base = bad_base[~bad_base['DB Spectrum'].isin(correct['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(correct['DB Spectrum'])]
    #Remove duplicates found in no_int
wrong_sub = wrong_sub[~wrong_sub['DB Spectrum'].isin(no_int['DB Spectrum'])]
bad_loc = bad_loc[~bad_loc['DB Spectrum'].isin(no_int['DB Spectrum'])]
ambiguous = ambiguous[~ambiguous['DB Spectrum'].isin(no_int['DB Spectrum'])]
other_mod = other_mod[~other_mod['DB Spectrum'].isin(no_int['DB Spectrum'])]
no_mod = no_mod[~no_mod['DB Spectrum'].isin(no_int['DB Spectrum'])]
bad_base = bad_base[~bad_base['DB Spectrum'].isin(no_int['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(no_int['DB Spectrum'])]
    #Remove duplicates found in wrong_sub
bad_loc = bad_loc[~bad_loc['DB Spectrum'].isin(wrong_sub['DB Spectrum'])]
ambiguous = ambiguous[~ambiguous['DB Spectrum'].isin(wrong_sub['DB Spectrum'])]
other_mod = other_mod[~other_mod['DB Spectrum'].isin(wrong_sub['DB Spectrum'])]
no_mod = no_mod[~no_mod['DB Spectrum'].isin(wrong_sub['DB Spectrum'])]
bad_base = bad_base[~bad_base['DB Spectrum'].isin(wrong_sub['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(wrong_sub['DB Spectrum'])]
    #Remove duplicates found in bad_loc
ambiguous = ambiguous[~ambiguous['DB Spectrum'].isin(bad_loc['DB Spectrum'])]
other_mod = other_mod[~other_mod['DB Spectrum'].isin(bad_loc['DB Spectrum'])]
no_mod = no_mod[~no_mod['DB Spectrum'].isin(bad_loc['DB Spectrum'])]
bad_base = bad_base[~bad_base['DB Spectrum'].isin(bad_loc['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(bad_loc['DB Spectrum'])]
    #Remove duplicates found in ambiguous
other_mod = other_mod[~other_mod['DB Spectrum'].isin(ambiguous['DB Spectrum'])]
no_mod = no_mod[~no_mod['DB Spectrum'].isin(ambiguous['DB Spectrum'])]
bad_base = bad_base[~bad_base['DB Spectrum'].isin(ambiguous['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(ambiguous['DB Spectrum'])]
    #Remove duplicates found in other_mod
no_mod = no_mod[~no_mod['DB Spectrum'].isin(other_mod['DB Spectrum'])]
bad_base = bad_base[~bad_base['DB Spectrum'].isin(other_mod['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(other_mod['DB Spectrum'])]
    #Remove duplicates found in no_mod
bad_base = bad_base[~bad_base['DB Spectrum'].isin(no_mod['DB Spectrum'])]
noid = noid[~noid['DB Spectrum'].isin(no_mod['DB Spectrum'])]
    #Remove duplicates found in bad_base
noid = noid[~noid['DB Spectrum'].isin(bad_base['DB Spectrum'])]


#%%
#Make the Sankey
n_noid = len(noid)
n_bad_base = len(bad_base)
n_no_mod = len(no_mod)
n_other_mod = len(other_mod)
n_ambiguous = len(ambiguous)
n_bad_loc = len(bad_loc)
n_wrong_sub = len(wrong_sub)
n_no_int = len(no_int)
n_correct = len(correct)
n_total = len(shared)

sankey_values =[n_correct,n_no_int,n_wrong_sub,n_bad_loc,
                      n_ambiguous,n_other_mod,n_no_mod,n_bad_base,n_noid,len(shared)]
sankey_destinations = [i for i in range(0,9)]
sankey_source = [9]*len(sankey_destinations)
labels = [f'Correct:{n_correct}',f'No Intensity:{n_no_int}',f'Wrong Substitution:{n_wrong_sub}',f'Bad Localization:{n_bad_loc}',
          f'Ambiguous Substitution:{n_ambiguous}',f'Other Modification:{n_other_mod}',f'No Modification:{n_no_mod}',
          f'Other Base:{n_bad_base}',f'No ID:{n_noid}','All']

fig = go.Figure(data=[go.Sankey(
    node=dict(
        label=labels,
    ),
    link=dict(
        source=sankey_source,
        target=sankey_destinations,
        value=sankey_values,
    )
)])
fig.update_layout(width=500,height=600,margin=dict(l=10, r=10, t=10, b=10))
fig.show()
fig.write_html(os.path.join(outputdir,"sankey_fates.html"))

#Figure with % labels
sankey_values_normalized = [x/n_total for x in sankey_values]
labels_percent = ['Correct:{:.1f}%'.format(sankey_values_normalized[0]*100),
'No Intensity:{:.1f}%'.format(sankey_values_normalized[1]*100),
'Wrong Substitution:{:.1f}%'.format(sankey_values_normalized[2]*100),
'Bad Localization:{:.1f}%'.format(sankey_values_normalized[3]*100),
'Ambiguous Substitution:{:.1f}%'.format(sankey_values_normalized[4]*100),
'Other Modification:{:.1f}%'.format(sankey_values_normalized[5]*100),
'No Modification:{:.1f}%'.format(sankey_values_normalized[6]*100),
'Other Base:{:.1f}%'.format(sankey_values_normalized[7]*100),
'No ID:{:.1f}%'.format(sankey_values_normalized[8]*100),
'All']

fig = go.Figure(data=[go.Sankey(
    node=dict(
        label=labels_percent,
    ),
    link=dict(
        source=sankey_source,
        target=sankey_destinations,
        value=sankey_values,
    )
)])
fig.update_layout(width=500,height=600,margin=dict(l=10, r=10, t=10, b=10))
fig.show()
fig.write_html(os.path.join(outputdir,"sankey_fates_percent.html"))


#%%
#Compare scores

shared['Score Shift'] = shared['DB Probability'] - shared['Probability']
plt.hist(shared['Score Shift'],bins=30)
plt.xlabel('Probability (DB - Discovery)')
plt.ylabel('Number of PSMs')
plt.title('All Target PSMs')
plt.savefig(os.path.join(outputdir,'Score Shift All PSMs.png'))

correct['Score Shift'] = correct['DB Probability'] - correct['Probability']
plt.figure()
plt.hist(correct['Score Shift'],bins=30)
plt.xlabel('Probability (DB - Discovery)')
plt.ylabel('Number of PSMs')
plt.title('Correctly identified PSMs')
plt.savefig(os.path.join(outputdir,'Score Shift Discovered PSMs.png'))


#%%
#Look at % ID by different variables
shared['Is Discovered'] = shared['Spectrum'].isin(correct['Spectrum'])

#By position
sharedbyPos = shared.groupby(by=['Substitution Peptide Position','Peptide Length'])['Is Discovered'].apply(lambda x: x.sum()/x.count())
widebyPos = sharedbyPos.reset_index().pivot_table(index='Substitution Peptide Position',
                                                  columns = 'Peptide Length',values='Is Discovered')
widebyPos = widebyPos.sort_index(ascending=False)
sns.heatmap(widebyPos,cmap = 'plasma_r',cbar=True, vmin=0, vmax=1)
plt.savefig(os.path.join(outputdir,'Discovered efficiency by position.png'))

#By type is weird
#Some of these PSMs could represent multiple substitution types (from different E. coli bases)
#We're just going to take the first option and run with that.
#Merge in substitution type info
def get_sub_type(pep):
    origin,destination = targetsbybase.get_group(pep)['Substitution Types'].iloc[0][0]   #The [0] takes the first of any types
    return origin,destination
shared[['Origin','Destination']] = shared['DB Peptide'].apply(lambda x: pd.Series(get_sub_type(x)))
shared['Origin'] = shared['Origin'].replace(p.dict_123ambiguous)
shared['Destination'] = shared['Destination'].replace(p.dict_123ambiguous)

sharedbyType = shared.groupby(by=['Origin','Destination'])['Is Discovered'].apply(lambda x: x.sum()/x.count())
widebyType = sharedbyType.reset_index().pivot_table(index='Origin', columns = 'Destination',values='Is Discovered')
widebyType = widebyType.sort_index(ascending=True)

plt.figure()
sns.heatmap(widebyType,cmap = 'plasma_r',cbar=True,vmin=0, vmax=1)
plt.savefig(os.path.join(outputdir,'Discovered efficiency by type.png'))


#%% 
#Look at intensity and dilution series stuff

#N targets by dilution
dil_n = shared.groupby(by=['Sample'])['Is Discovered'].sum()
dil_n = dil_n.sort_index(ascending=False)
plt.bar(dil_n.index.astype(str),dil_n)
plt.ylabel('# of Target Spectra')
plt.xlabel('% S. typhi')
plt.savefig(os.path.join(outputdir,'Bar N Targets.png'))


#Delta Quant
shared['Intensity shift'] = shared['DB Intensity'] - shared['Intensity']
mask_intense_sub = shared['Intensity'] > 0
mask_intense_db = shared['DB Intensity'] > 0
mask_intense_both = mask_intense_db & mask_intense_sub
shared['log(Intensity shift)'] = np.log10(shared['Intensity shift'].replace(0,np.nan))
shared.loc[shared['log(Intensity shift)']<0,'log(Intensity shift)'] = 0 #Pseudolog
shared.loc[~mask_intense_db,'log(Intensity shift)'] = (-1)*np.log10(shared.loc[~mask_intense_db,'Intensity shift'].replace(0,np.nan)*(-1))

plt.figure()
plt.hist(shared.loc[mask_intense_both,'log(Intensity shift)'],bins=20,histtype='step',label='Shifted')
plt.hist(shared.loc[~mask_intense_sub,'log(Intensity shift)'],bins=20,histtype='step',label='Unquantified Sub')
plt.hist(shared.loc[~mask_intense_db,'log(Intensity shift)'],bins=20,histtype='step',label='Unquantified DB')
plt.bar(x=0,height=(shared['Intensity shift'] == 0).sum(),width=0.2,color='k',label='No shift')
plt.xlabel('Intensity shift (DB-discovery) psuedolog')
plt.ylabel('# of Target Spectra')
plt.legend()
plt.savefig(os.path.join(outputdir,'Intensity shift.png'))


#Dynamic range
range_undisc = shared[~shared['Is Discovered']].groupby(by=['DB Peptide'])['Intensity'].apply(lambda x: x.max()-x.replace(0,np.nan).min())
range_disc = shared[shared['Is Discovered']].groupby(by=['DB Peptide'])['Intensity'].apply(lambda x: x.max()-x.replace(0,np.nan).min())
undisc_0 = range_undisc == 0
disc_0 = range_disc == 0
range_undisc = np.log10(range_undisc.replace(0,np.nan))
range_undisc.loc[undisc_0] = 0
range_disc = np.log10(range_disc.replace(0,np.nan))
range_disc.loc[disc_0] = 0

plt.figure()
plt.hist(range_disc,bins=20,histtype='step',label='Discovered Targets')
plt.hist(range_undisc,bins=20,histtype='step',label='Undiscovered Targets')
plt.xlabel('PSM Dynamic Range (max-min) pseudolog')
plt.ylabel('# of Target Peptides')
plt.legend()
plt.savefig(os.path.join(outputdir,'Target peptide dynamic range.png'))

#%%
def dilution_dropout(row):
    """Steps through a descending dilution series
    Allows 1 'strike' of missing value
    Returns the Sample/column lowest in the series with an intensity value"""

    last_non_na_col = None
    cumulative_na_count = 0
    # Iterate over each element in the row
    for col, val in row.items():
        # Check if the value is not NaN
        if not pd.isna(val):
            # Update the last_non_na_col with the current column
            last_non_na_col = col
        else:
            # Increment the consecutive NaN count
            cumulative_na_count += 1
            if cumulative_na_count == 2:        #If this is the second strike
                return last_non_na_col          #Return the previous col with a value
            
    return last_non_na_col                      #If we make it through whole series without a second strike, return last col

#Last observed quantification
dil_disc = shared[shared['Is Discovered']].groupby(by=['Sample','DB Peptide'])['Intensity'].apply(lambda x: x.replace(0,np.nan).mean())
dil_undisc = shared[~shared['Is Discovered']].groupby(by=['Sample','DB Peptide'])['Intensity'].apply(lambda x: x.replace(0,np.nan).mean())
dil_disc = np.log10(dil_disc.replace(0,np.nan))
dil_undisc = np.log10(dil_undisc.replace(0,np.nan))
dil_disc = dil_disc.reset_index().pivot_table(columns='Sample',index='DB Peptide',values='Intensity')
dil_undisc = dil_undisc.reset_index().pivot_table(columns='Sample',index='DB Peptide',values='Intensity')
dil_disc.sort_index(ascending=False,axis=1,inplace=True)
dil_undisc.sort_index(ascending=False,axis=1,inplace=True)

dil_disc['Last Quantified Sample'] = dil_disc.apply(dilution_dropout,axis=1)
dil_undisc['Last Quantified Sample'] = dil_undisc.apply(dilution_dropout,axis=1)

distributions = []
for sample in samples:
    try:
        distributions += [dil_disc.loc[dil_disc['Last Quantified Sample'] == sample,sample].to_list()]
    except KeyError:
        distributions += [7]
fig,ax = plt.subplots()
plt.boxplot(distributions)
plt.ylabel('log Intensity')
ax.set_xticklabels(order)
plt.xlabel('% S. typhi')
plt.savefig(os.path.join(outputdir,'Last Seen Intensity Distribution.png'))

#%%
#Export

shared.to_csv(os.path.join(outputdir,'Shared DB_Discovery Table with Filters.csv'),index=False)
noid.to_csv(os.path.join(outputdir,'NoID in discovery.csv'),index=False)
bad_base.to_csv(os.path.join(outputdir,'Base base in discovery.csv'),index=False)
no_mod.to_csv(os.path.join(outputdir,'No modification in discovery.csv'),index=False)
other_mod.to_csv(os.path.join(outputdir,'Other modification in discovery.csv'),index=False)
ambiguous.to_csv(os.path.join(outputdir,'Ambiguous with PTM in discovery.csv'),index=False)
bad_loc.to_csv(os.path.join(outputdir,'Bad localization in discovery.csv'),index=False)
wrong_sub.to_csv(os.path.join(outputdir,'Wrong substitution in discovery.csv'),index=False)
no_int.to_csv(os.path.join(outputdir,'No intensity in discovery.csv'),index=False)
correct.to_csv(os.path.join(outputdir,'Correct in discovery.csv'),index=False)
p.timestamp('Done')
import winsound
duration = 750  # milliseconds
freq = 500  # Hz
winsound.Beep(freq, duration)
winsound.Beep(freq, duration)

# %%
