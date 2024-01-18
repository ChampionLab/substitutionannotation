#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:48:37 2021
Note: This script had hardcoded paths. These must be replaced, and all are at the top:
OutputDirectory, and two input TSVs of tryptic peptide lists imported into the variables dfsalty and dfecoli.
@author: taylorlundgren
"""

"""Search two peptide lists and identify singly substituted peptides (SSP) between them.
This will EXCLUDE 
    ...exactly similar peptide sequences
    ...sequences that match but have an early truncation, with or without SSP
    ...substitutions without K/R as the last amino acid (new cutsite means it is not SSP)
This will INCLUDE
    ...substitutions of K/R NOT the last sequence (these represent novel SSP with one missed cleavage)
"""
import pandas as pd
from tqdm import tqdm
import os

#%%
#Import tsv list of tryptic peptides from 2 organisms. Note: include 1 missed cleavage to identify single scissor substitutions that remove a cutsite.
    #Include 2 missed cleavage peptides to identify doubly substitued peptides that remove two cutsites. 
dfsalty = pd.read_csv(r"PathTo\SALTY_TrypticPeptides.tsv", sep='\t',usecols=['Base Sequence'])
dfecoli = pd.read_csv(r"PathTo\Ecoli_tryptic_peptides.tsv", sep='\t',usecols=['Base Sequence'])
dfecoli = dfecoli.drop_duplicates(subset='Base Sequence')
dfsalty = dfsalty.drop_duplicates(subset='Base Sequence')

OutputDirectory = ""

#%%
def find_SSP(str1, str2):
    """
    A function to see if str2 is a SSP of str1.
    'Normal' substitutions have the same length and one difference.
    Scissor substitutions can only be found by including peptides with 1 missed cleavage in the comparison.
    Scissor substitutions will be identified by comparing a peptide with 0 missed cleavages to another with 1 missed cleavage.
    """
    if len(str1 != len(str2)):  #Keep comparisons to tryptic peptides of same length
        return None

    diff_indices = [i for i, (c1, c2) in enumerate(zip(str1, str2)) if c1 != c2]
    if len(diff_indices) == 1:                              #Check for SINGLY substituted sequence
        if len(str1) == len(str2):                          #Logic for same length peptide sequences
            if diff_indices[0]+1 == len(str1):                  #If peptide C-term aa is different
                if str2[diff_indices[0]] in ['K','R']:
                    return str2                                 #Return SSP if tryptic cleavage rule is preserved
                else:
                    return None                                 #Prevent matching SSP to semi-enzymatic peptide fragment
            else:
                return str2                                 #Return sequence representing SSP

    else:
        return None

#%%
def find_DSP(str1, str2):
    """
    A function to see if str2 is a DSP of str1.
    'Normal' substitutions have the same length and one difference.
    Single scissor substitutions can only be found by including peptides with 1 missed cleavage in the comparison.
    Double scissor substitutions can only be found by including peptides with 2 missed cleavage in the comparison.
    """
    if len(str1 != len(str2)):  #Keep comparisons to tryptic peptides of same length
        return None

    diff_indices = [i for i, (c1, c2) in enumerate(zip(str1, str2)) if c1 != c2]
    if len(diff_indices) == 2:                              #Check for DOUBLY substituted sequence
        if len(str1) == len(str2):                          #Logic for same length peptide sequences
            if diff_indices[1]+1 == len(str1):                  #If peptide C-term aa is different
                if str2[diff_indices[1]] in ['K','R']:
                    return str2                                 #Return SSP if tryptic cleavage rule is preserved
                else:
                    return None                                 #Prevent matching SSP to semi-enzymatic peptide fragment
            else:
                return str2                                 #Return sequence representing SSP

    else:
        return None



#%%
#Find salty peptides representing a ssp from ecoli peptides
print('Finding SSP from first list')
progress_bar = tqdm(total=len(dfecoli))
dfecoli['SALTY SSP Sequence'] = dfecoli['Base Sequence'].apply(lambda x: (list((find_SSP(x, y) for y in dfsalty['Base Sequence'] if find_SSP(x, y))),progress_bar.update(1))[0])
dfecoli['Is SSP'] = dfecoli['SALTY SSP Sequence'].astype(bool)
progress_bar.close()
#%%
#Find ecoli peptides representing a ssp from salty peptides
print('Finding SSP from second list')
progress_bar = tqdm(total=len(dfsalty))
dfsalty['ECOLI SSP Sequence'] = dfsalty['Base Sequence'].apply(lambda x: (list((find_SSP(x, y) for y in dfecoli['Base Sequence'] if find_SSP(x, y))),progress_bar.update(1))[0])
dfsalty['Is SSP'] = dfsalty['ECOLI SSP Sequence'].astype(bool)
progress_bar.close()
#%%
#Find salty peptides representing a dsp from ecoli peptides
print('Finding DSP from first list')
progress_bar = tqdm(total=len(dfecoli))
dfecoli['SALTY DSP Sequence'] = dfecoli['Base Sequence'].apply(lambda x: (list((find_DSP(x, y) for y in dfsalty['Base Sequence'] if find_DSP(x, y))),progress_bar.update(1))[0])
dfecoli['Is DSP'] = dfecoli['SALTY DSP Sequence'].astype(bool)
progress_bar.close()
#%%
#Find ecoli peptides representing a dsp from salty peptides
print('Finding DSP from second list')
progress_bar = tqdm(total=len(dfsalty))
dfsalty['ECOLI DSP Sequence'] = dfsalty['Base Sequence'].apply(lambda x: (list((find_DSP(x, y) for y in dfecoli['Base Sequence'] if find_DSP(x, y))),progress_bar.update(1))[0])
dfsalty['Is DSP'] = dfsalty['ECOLI DSP Sequence'].astype(bool)
progress_bar.close()

#%%
dfecoli['Is Exact'] = dfecoli['Base Sequence'].isin(dfsalty['Base Sequence'])
dfsalty['Is Exact'] = dfsalty['Base Sequence'].isin(dfecoli['Base Sequence'])

#Report the fraction of SSPs
dffrac = pd.DataFrame(index=['From ECOLI','From SALTY'])
dffrac.loc['From SALTY','Total'] = len(dfsalty)
dffrac.loc['From SALTY','Single Sub'] = dfsalty['Is SSP'].sum()
dffrac.loc['From SALTY','Double Sub'] = dfsalty['Is DSP'].sum()
dffrac.loc['From SALTY','Single AND Double'] = (dfsalty['Is DSP'] & dfsalty['Is SSP']).sum()
dffrac.loc['From SALTY','Single AND Exact'] = (dfsalty['Is SSP'] & dfsalty['Is Exact']).sum()
dffrac.loc['From SALTY','Double AND Exact'] = (dfsalty['Is DSP'] & dfsalty['Is Exact']).sum()
dffrac.loc['From SALTY','Single AND Double AND Exact'] = (dfsalty['Is DSP'] & dfsalty['Is SSP'] & dfsalty['Is Exact']).sum()
dffrac.loc['From SALTY','Exact Match'] = dfsalty['Is Exact'].sum()
dffrac.loc['From SALTY','Different Sequences'] = (dffrac.loc['From SALTY','Total'] - 
                                                dffrac.loc['From SALTY','Single Sub']
                                                - dffrac.loc['From SALTY','Double Sub'] + dffrac.loc['From SALTY','Single AND Double']
                                                - dffrac.loc['From SALTY','Exact Match'] + dffrac.loc['From SALTY','Single AND Exact'] 
                                                    + dffrac.loc['From SALTY','Double AND Exact'] - dffrac.loc['From SALTY','Single AND Double AND Exact'])

dffrac.loc['From ECOLI','Total'] = len(dfecoli)
dffrac.loc['From ECOLI','Single Sub'] = dfecoli['Is SSP'].sum()
dffrac.loc['From ECOLI','Double Sub'] = dfecoli['Is DSP'].sum()
dffrac.loc['From ECOLI','Exact Match'] = dfecoli['Is Exact'].sum()
dffrac.loc['From ECOLI','Single AND Double'] = (dfecoli['Is DSP'] & dfecoli['Is SSP']).sum()
dffrac.loc['From ECOLI','Single AND Exact'] = (dfecoli['Is SSP'] & dfecoli['Is Exact']).sum()
dffrac.loc['From ECOLI','Double AND Exact'] = (dfecoli['Is DSP'] & dfecoli['Is Exact']).sum()
dffrac.loc['From ECOLI','Single AND Double AND Exact'] = (dfecoli['Is DSP'] & dfecoli['Is SSP'] & dfecoli['Is Exact']).sum()
dffrac.loc['From ECOLI','Exact Match'] = dfecoli['Is Exact'].sum()
dffrac.loc['From ECOLI','Different Sequences'] = (dffrac.loc['From ECOLI','Total'] - 
                                                dffrac.loc['From ECOLI','Single Sub']
                                                - dffrac.loc['From ECOLI','Double Sub'] + dffrac.loc['From ECOLI','Single AND Double']
                                                - dffrac.loc['From ECOLI','Exact Match'] + dffrac.loc['From ECOLI','Single AND Exact'] 
                                                    + dffrac.loc['From ECOLI','Double AND Exact'] - dffrac.loc['From ECOLI','Single AND Double AND Exact'])

#%%
#Export tables
dffrac.to_csv(os.path.join(OutputDirectory,'Ecoli_v_Salty_Similarity_Summary.csv'))
dfsalty.to_csv(os.path.join(OutputDirectory,'All_SALTY_to_ECOLI.csv'))
dfecoli.to_csv(os.path.join(OutputDirectory,'All_ECOLI_to_SALTY.csv'))

dff = dfsalty[dfsalty['Is SSP']]
dff.rename({'Base Sequence':'SALTY Sequence'})
dff.to_csv(os.path.join(OutputDirectory,'SSP_SALTY_to_ECOLI.csv'))

dff = dfecoli[dfecoli['Is SSP']]
dff.rename({'Base Sequence':'ECOLI Sequence'})
dff.to_csv(os.path.join(OutputDirectory,'SSP_ECOLI_to_SALTY.csv'))

dff = dfsalty[dfsalty['Is DSP']]
dff.rename({'Base Sequence':'SALTY Sequence'})
dff.to_csv(os.path.join(OutputDirectory,'DSP_SALTY_to_ECOLI.csv'))

dff = dfecoli[dfecoli['Is DSP']]
dff.rename({'Base Sequence':'ECOLI Sequence'})
dff.to_csv(os.path.join(OutputDirectory,'DSP_ECOLI_to_SALTY.csv'))

# %%
