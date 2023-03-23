#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:48:37 2021

@author: taylorlundgren
"""

"""Search two peptide lists and identify singly substituted peptides (SSP) between them.
This will EXCLUDE 
    ...exactly similar peptide sequences
    ...peptides which have an early truncation, with out without SSP
    ...substitutions of K/R as the last sequence (new cutsite means it is not SSP)
This will INCLUDE
    ...substitutions of K/R NOT the last sequence (these represent novel SSP with one missed cleavage)
"""
import pandas as pd
import numpy as np
import re

def replaceChar(string,position):  #Replaces the character at index position with a regex wildcard
    result = string[: position] + '.' + string[position + 1:]
    return result

def GetSubstitution(indx,originAA,destinationpep):
    destinationAA = destinationpep[i]
    result = [originAA+'->'+destinationAA]
    return result

#Import all peptides in each database, set up boolean vector to determine SSPs
dfsalty = pd.read_csv('/Users/taylorlundgren/PythonLab/Salty_Peptides_1.tsv', sep='\t',usecols=['Base Sequence'])
dfecoli = pd.read_csv('/Users/taylorlundgren/PythonLab/Coli_Peptides_1.tsv', sep='\t',usecols=['Base Sequence'])
dfecoli['IsSSP']=[False]*len(dfecoli)
dfsalty['IsSSP'] = [False]*len(dfsalty)

#Populate empty list column for listing SSP sequence in other DB, represented substitutions
dfsalty['EColi SSP Peptide'] = np.empty((len(dfsalty), 0)).tolist()
dfsalty['Substitution Types'] = np.empty((len(dfsalty), 0)).tolist()


i=0
for peptide in dfsalty['Base Sequence']:
    if peptide[-1] != ('R' or 'K'):     #exclude substitutions involving K/R at the C terminus
        indx = dfsalty.index[dfsalty['Base Sequence'] == peptide]
        
        if (dfecoli['Base Sequence']==peptide).any():   #Set exact matches to False
            dfsalty.loc[indx,'IsSSP'] = False
            dfsalty.loc[indx,'EColi SSP Peptide'] = []
            dfsalty.loc[indx,'Substitution Types'] = []
            continue                                    #And skip looking for SSPs
            
        for i in range(0,len(peptide)):     #Iterate for each position in the peptide sequence
            regexstring = replaceChar(peptide,i)    #Float one character to match to any other character
                            
            #if the full match is true, or if there is already an identified SSP, set value to true
            dfsalty.loc[indx,'IsSSP'] = dfsalty.loc[indx,'IsSSP'] | bool(dfecoli['Base Sequence'].apply(lambda x: re.fullmatch(regexstring,x)).any())
            #Get the corresponding peptides which were ID'd as SSPs
            dfsalty.loc[indx,'EColi SSP Peptide'] = dfsalty.loc[indx,'EColi SSP Peptide'] + list(dfecoli['Base Sequence'].apply(lambda x: re.fullmatch(regexstring,x)).any())
            #Get the corresponding substitution types
            dfsalty.loc[indx,'Substitution Types'] = dfsalty.loc[indx,'Substitution Types'] + (
                [GetSubstitution(i,peptide[i],item) for item in list(dfecoli['Base Sequence'].apply(lambda x: re.fullmatch(regexstring,x)).any())]
                )

        
            
#Export the portion containing only SSPs

dfSSP = dfsalty[dfsalty['IsSSP']]
dfSSP.rename({'Base Sequence':'Salty Sequence'})
dfsalty.to_csv('/Users/taylorlundgren/opt/anaconda3/lib/python3.8/tjl/params/SALTYagainstECOLIsspUnfiltered.csv')
dfSSP.to_csv('/Users/taylorlundgren/opt/anaconda3/lib/python3.8/tjl/params/SALTYagaisntECOLIssp.csv')