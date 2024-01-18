#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:43:35 2021
Note: this code has hardcoded paths for inputdir, dfdm, dfdangermods
@author: taylorlundgren
"""
"""Program to count the subs from FragPipe v17,
 Mass-Offset-Substitutions Good LFQ workflow output"""
#Import packages used
import pandas as pd
import numpy as np
import os
import re 
import warnings
import datetime
#%%
#~~~~~These are the things that may need to be changed run by run~~~~~~~~~~~~~~~~~

#Change the input directory to the output folder of FragPipe
inputdir = r'C:\Users\taylo\Documents\PythonLab\MixedExperiment\Final\eColi'

#Summaries by sample extract the sample name from the full file name in psm.tsv
#The regex to do so changes based on the LC system used. Change the match pattern
#below to properly extract the sample information from the full file name
#Be sure to prefix the string with an r -- i.e. r'Regex match pattern'
regSampleMatch = r'(.*)_Slot'   #Works with MASSIVE ######
# Example I've used in my evosep samples = r'(.*(?<=_S\d-).*)_'


#%%
#Reference items

outputdir = os.path.join(inputdir,'FindSubsOutput')
dfdm = pd.read_csv(r'C:\Users\taylo\anaconda3\Lib\tjl\params\dmMatrix.csv',index_col=0)
dfdangermods = pd.read_csv(r'C:\Users\taylo\anaconda3\Lib\tjl\params\dangermods.csv')

#Columns to import from Fragger's psm.tsv output
columns = ['Spectrum','Peptide','Modified Peptide','Calibrated Observed Mass',
           'Intensity','Peptide Length','Assigned Modifications','Protein',
           'Hyperscore','PeptideProphet Probability','Retention',
           'Ion Mobility','Calibrated Observed M/Z','Number of Missed Cleavages']

#Use a dictionary which assigns I and L to 'Leu/Ile' to reflect isobaric masses
dict_123 = {'G': 'Gly',
 'A': 'Ala',
 'S': 'Ser',
 'P': 'Pro',
 'V': 'Val',
 'T': 'Thr',
 'L': 'Leu/Ile',
 'I': 'Leu/Ile',
 'N': 'Asn',
 'D': 'Asp',
 'Q': 'Gln',
 'K': 'Lys',
 'E': 'Glu',
 'M': 'Met',
 'H': 'His',
 'F': 'Phe',
 'R': 'Arg',
 'C': 'Cys',
 'Y': 'Tyr',
 'W': 'Trp'}

#Dictionary for 3 letter code to one letter code
dict_321 = {"Gly": 'G', 
            "Ala" : 'A', 
            "Ser" : 'S', 
            "Pro" : 'P', 
            "Val" : 'V', 
            "Thr" : 'T', 
            "Ile" : 'I',
            "Leu" : "L",
            "Asn" : 'N', 
            "Asp" : 'D', 
            "Gln" : 'Q', 
            "Lys" : 'K', 
            "Glu" : 'E', 
            "Met" : 'M', 
            "His" : 'H',
            "Phe" : 'F', 
            "Arg" : 'R', 
            "Cys" : 'C', 
            "Tyr" : 'Y',
            "Trp" : 'W',
            "Cis" : 'C'
            }

#%%
def getPTMs(row):
    modAAs = re.findall('([A-Z])\(',row['Assigned Modifications'])
    dms= re.findall('\((.*?)\)',row['Assigned Modifications'])
    modpos= re.findall('(\d+)[A-Z]\(',row['Assigned Modifications'])
    #Set tolerance to 25ppm -- this should be ~2x the expected mass deviation,
    #as we extrapolate from the assigned delta mass which was within tolerance
    #and do not have the 'true' delta mass
    atol = row['Calibrated Observed Mass'] *2.5*10**(-5)
    ptms=[]
    if modAAs:
        for res,dm,position in zip(modAAs,dms,modpos):
            safebool = True
            res3 = dict_123.get(res)
            if res3 is None:
                warnings.warn('Could not find 3 letter code for '+res)
                return None
            #Get index for SAAV which match the observed delta mass at the observed residue    
            subindx = np.isclose(dfdm[res3].astype(float),float(dm),atol=atol,rtol=0) 
            
            #Get the index of common modifications which match the delta mass
            dangerindx = np.isclose(dfdangermods['Mass Shift'],float(dm),atol=atol,rtol=0)
            
            #Apply conditional filters for common modifications
            #Check if modified residue is explicitly listed in the identified danger mods
            if dfdangermods.loc[dangerindx,'Modified Residues'].str.contains(res).any():
                safebool = False
            #If the modified residue doesn't explicitly apply to possible danger mods, check ambiguity
            if dfdangermods.loc[dangerindx,'Modified Residues'].str.contains(res).any()==False:
                if dfdangermods.loc[dangerindx,'Modified Residues'].str.contains('B').any():
                    safebool = False #B is used for mods which apply to any residue
                if dfdangermods.loc[dangerindx,'Modified Position'].any():#If the mod applies to a position
                    for possiblepos in dfdangermods.loc[dangerindx,'Modified Position']:
                        if possiblepos == float(position):
                            safebool = False
            
            #Clear dangermod index if it failed to meet aa identify/position criteria
            if safebool: dangerindx = len(dangerindx)*[False]
               
            
            if (dfdm[subindx][res3].any() & bool(dfdangermods[dangerindx]['Modification'].any())):
                ptms += [res3+'->'+dfdm[subindx][res3].index[0]+' or '+dfdangermods.loc[dangerindx,'Modification'].values[0]]
            elif dfdm[subindx][res3].any():
                ptms += [res3+'->'+dfdm[subindx][res3].index[0]]
            elif dfdangermods[dangerindx]['Modification'].any():
                ptms += [dfdangermods.loc[dangerindx,'Modification'].values[0]]
    if ptms:
        return ','.join(ptms)
    else:
        return None
#%%
def GetBaseBy(spec,pep,byboth,by):
    """

    Parameters
    ----------
    spec : str
        Filename.
    pep : str
        Peptide sequence.
    byboth : pandas DF
        Value grouped by peptide sequence AND filename.
    by : pandas DF
        value grouped by sequence only.

    Returns
    -------
        value of genomic cognate from the same file if available, from any file
        otherwise, and None if not found in any file.

    """
    indx = (spec,pep)
    try:
        result = byboth.loc[indx]      #Try to get intensity matching file AND peptide
    except:
       # print('Failed to get groupby result for '+spec+' and '+pep)
        try: 
            result = by.loc[pep]    #Try to get intensity matching only peptide
        except:
            return None                    #No match returns none
    if result == 0:                        #If first attempt found 0 intensity, try just peptide match
        result = by.loc[pep]
    if result == 0:                        #If even matching only by peptide returns 0 intensity, return none
            return None
    return result
#%%
def GetBaseTrue(spec,pep,byboth):
    """
    Takes a file name, peptide sequence, and series whose index is a 
    tuple of file name, peptide sequence.
    Returns True if the injection/peptide index corresponds to an intensity
    greater than 0, False otherwise.
    """
    indx = (spec,pep)
    try: 
        result = byboth.loc[indx]
    except:
        return False
    if result > 0:
        return True   
    else:
        return False

#%%
def FindModSequence(dfin):
    """Given a peptide, mods and possible substitution sites, output possible sequences
    found
    Returns a dataframe with the Origin (genomic cognate) aa,
    Destination (substituted variant) aa,
    and Modified Sequences (entire peptide with substituted sequence,
                            two sequences if I/L ambiguous)
    """
    df = dfin.copy()
    #Match origin and destination from PTMs column from getPTMs 
    df['Origin'] = df['PTMs'].str.extract('(.{3})(?=->)')
    df['Destination'] = df['PTMs'].str.extract('(?<=->)(.{3})')
    #Location of I/L ambiguous substitutions
    maskIL = df['Destination'] == 'Leu'
    
    #Make new sequence string
    df['Modified Sequences'] = (df['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        df['Destination'].replace(dict_321) +
                        df['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0])
    #Subset I/L ambiguous sequences
    dff = df[maskIL]
    #Make list of two possible sequences
    dff.loc[:,'Modified Sequences'] = (dff['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        'I' +
                        dff['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0] + ','+
                        dff['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        'L' +
                        dff['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0])
    #Write ambiguous sequences over single outputs
    df.loc[maskIL,'Modified Sequences'] = dff['Modified Sequences']
                        
    return df[['Origin','Destination','Modified Sequences']]    
    
    
#%%
def ExpIL(dfin,sequencecolumn): 
    """"Split list of possible sequences in provided
 df[sequencecolumn using , as delimiter, expanding the rows for each sequence."""
    df = dfin.copy()
    df[sequencecolumn] = df.loc[:,sequencecolumn].str.split(',')
    df = df.explode(sequencecolumn)
    return df

#%%
#Import data
df = pd.read_csv(os.path.join(inputdir,'combined_protein.tsv'), sep = '\t')

if os.path.isdir(outputdir):
    print('Output directory already exists -- code will overwrite existing files with common names')
else: 
    os.mkdir(outputdir)
    print('Created output directory '+outputdir)


#Get filenames

filenames = pd.read_csv(os.path.join(inputdir,'filelist_ionquant.txt'),sep='\t')
#Drop spectra directory and flag column
filenames = filenames.loc[~filenames['flag'].str.contains('specdir'),'value']
#Drop psm.tsv
filenames = filenames.replace(regex=r'(\\.*)',value='')
#Drop duplicated names
filenames.drop_duplicates(inplace=True)
filenames.to_csv(os.path.join(outputdir,'filenames.csv'),index=False)


#Pull out PSM info
print('Retreiving PSM info...')
now = datetime.datetime.now()
print(str(now))

df = pd.DataFrame()
#Import the psm.tsv of each sample
for file in filenames:
    dfiter = pd.read_csv(os.path.join(inputdir,file,'psm.tsv'), sep='\t', usecols=columns)
    df = pd.concat([df,dfiter])
df.reset_index(drop=True,inplace=True)    

#Sometimes fragger hangs up on filenames with a . and renames them to _
#This creates a discrepancy between the preserved filename under 'Spectrum' 
#in psm.tsv and the modified filenames used elsewhere by fragger
#This is the case for data in MASSIVE #####, so I reconcile the name here
df['Spectrum'] = df['Spectrum'].str.replace(r'SALTY3.5_','SALTY3_5_')

#Replace NaN values with 'None' for compatibility with labeling modifications
df['Assigned Modifications'].fillna('None',inplace=True)

#Align filename from fragpipe manifest, spectrum name from psm.tsv, 
#and sample name from fragpipe manifest
df['Full Spectrum'] = df['Spectrum']
df['Sample'] = df['Spectrum'].str.extract(regSampleMatch)


#Find assigned modifications which match a sub
print('Naming modifications, this may take some time...')
now = datetime.datetime.now()
print(str(now))
df['PTMs'] = df.apply(lambda x:getPTMs(x),axis=1)
df.loc[:,'PTMs'] = df['PTMs'].fillna('None')

#Extract decoy PSMs, if exported
print('Extracting decoy PSMs...')
now = datetime.datetime.now()
print(str(now))
dfdecoy = df.loc[df['Protein'].str.contains('rev_')]
df = df.loc[~df['Protein'].str.contains('rev_')]
dfdecoy.to_csv(os.path.join(outputdir,'DiscoveryDecoyPSMs.csv'),index =False)
#%%
#Report if a psm is a possible substitution 'Is Sub'
#or if a psm could have a different modification 'Is Danger'
print('Making binary filters...')
now = datetime.datetime.now()
print(str(now))

df['Is Sub'] = df['PTMs'].str.contains('->')
df['Is Danger'] = df['PTMs'].str.contains('or')

#Filter out 0 intensity PSM
#'Has Intensity' column for PSMs - True values have an intensity > 0, False do not.
df.loc[:,'Has Intensity'] = df['Intensity'].apply(bool)

# 'Is Mod' column is True for PSMs of any modification
df.loc[:,'Is Mod'] = ~df['PTMs'].str.contains('None')

#Find genomic cognate peptides to identified substitutions
#List of genomic cognate peptide sequences
modpeplist = df.loc[df['Is Sub'],'Peptide'].drop_duplicates()
#Column 'Is Base' is True if the psm sequence is a cognate AND the psm is unmodified
df.loc[:,'Is Base'] = df['Peptide'].isin(modpeplist) & ~df['Is Mod']

#Columns with a True value for psms that represent a confident substitution
#'Filtered' for substituted PSMs with intensity >0
#'Filtered (No Intensity) for substitued PSMs with intensity >=0
df.loc[:,'Filtered'] = ~df['Is Danger']&df['Has Intensity']&df['Is Sub']
df.loc[:,'Filtered (No Intensity)'] = ~df['Is Danger'] & df['Is Sub']
#A column full of True values, an admittedly dumb way of counting the number of PSMs but it works
df.loc[:,'Is Spectrum'] = [True]*len(df)
#%%
#Get Modified Sequence

df[['Origin','Destination','Modified Sequence']] = FindModSequence(df)
df = ExpIL(df,'Modified Sequence')
df['Modified Sequence'].fillna(df['Modified Peptide'],inplace=True)
df['Modified Sequence'].fillna(df['Peptide'],inplace=True)
#%%
#Append the genoic cognate Intensity, Retention, and Ion mobility information
#as columns next to each substituted PSM
print('Getting Base Peptide Info...')
now = datetime.datetime.now()
print(str(now))

#Genomic cognate information to add
characteristics = ['Intensity','Retention','Ion Mobility']
#Parse out substitutions
dfsubs = df[df['Is Sub']]

#Get genomic cognate information one at a time
for char in characteristics:
    #Get cognate information from within the same sample
    bybaseboth = df[df['Is Base']].groupby(by=['Spectrum','Peptide'])[char].max()
    #Get cognate information from any sample
    bybase = df[df['Is Base']].groupby(by=['Peptide'])[char].max()
    #Apply a function which returns cognate information first from the same sample,
    #then from any sample
    dfsubs[('Base '+char)] = dfsubs[['Spectrum','Peptide']].apply(
        lambda x:GetBaseBy(*x,bybaseboth,bybase),axis=1)
    if char == 'Intensity':
        #Report if the base information was from the same sample (True)
        #or from any sample (False)
        #in 'True Base' column
        dfsubs['True Base'] = dfsubs[['Spectrum','Peptide']].apply(
            lambda x:GetBaseTrue(*x,bybaseboth),axis=1)

#%%
#Output files
print('Generating Reports...')
now = datetime.datetime.now()
print(str(now))

#List of substituted peptide PSMs with or without intensity
dffilteredpsm = dfsubs[dfsubs['Filtered (No Intensity)']]
dffilteredpsm.to_csv(os.path.join(outputdir,'Filtered SSP PSM No Intensity Filter.csv'),index=False)
dffilteredpsm = dfsubs[dfsubs['Filtered']]
dffilteredpsm.to_csv(os.path.join(outputdir,'Filtered SSP PSM.csv'),index=False)

#Summary data output
metrics = ['Is Spectrum','Filtered','Filtered (No Intensity)','Is Danger','Has Intensity']
#Dictionary to rename columns to more verbose/readable outputs
dict_metric = { 'Is Spectrum':'Total # of PSMs',
               'Filtered':'Filtered SSP PSMs',
               'Filtered (No Intensity)':'Filtered SSP PSMs (No Intensity Filter)',
               'Is Danger':'# of Dangermod Ambiguous Sub PSMs',
               'Has Intensity':'# of Intense PSMs'}

#Get sum of each metric by sample. Make readable and export
dfsummary = df.groupby(by=['Sample'])[metrics].sum().T
dfsummary['Total'] = dfsummary.sum(axis=1)
dfsummary = dfsummary.rename(index=dict_metric)
dfsummary.to_csv(os.path.join(outputdir,'PSM Summary.csv'))

#Export all PSMs with added columns
df.to_csv(os.path.join(outputdir,'AllPSMsAndFilters.csv'),index=False)

#Get summary data by mod and by file
print('Parsing data by mod/file...')
now = datetime.datetime.now()
print(str(now))

#Each modification combination is reported
byfilemod = df.groupby(by=['Sample','PTMs'])['Is Spectrum'].count()
#Reshape table
dfbymod = byfilemod.reset_index().pivot(columns='Sample',index='PTMs',values='Is Spectrum')
dfbymod.loc[:,'Total']=dfbymod.sum(axis=1) 
dfbymod.to_csv(os.path.join(outputdir,'PSMs By Mod and File.csv'))  


print('Done.')
now = datetime.datetime.now()
print(str(now))



