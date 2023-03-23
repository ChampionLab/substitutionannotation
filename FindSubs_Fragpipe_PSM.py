#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 11:43:35 2021

@author: taylorlundgren
"""
"""Program to count the subs from FragPipe v17,
 Mass-Offset-Substitutions Good LFQ workflow output"""

import pandas as pd
import numpy as np
import os
import re 
import datetime
from tjl.params import tjlassets as p

#Regex capture pattern to get sample name from file name
#regexSample = 'GB_(.*)-'
    

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
#Reset Variables
#%reset -f
#%%
def getModaa(row):
    """
    Parameters
    ----------
    row : pandas Series, representing a row from a dataframe
    Contains Modified Peptide
    
    Returns 
    -------
    A series with a string indicating the smallest common denominator of possibly
    modified aa, and another string with every potentially modified aa
    """
    
    if '[' in(row['Modified Peptide']):
        
        unique_aa = [aa for row in row['MSFragger Localization']
                              for aa in row if aa.islower()]
        unique_aa = [letter.upper() for letter in unique_aa]
        return unique_aa
    else:
        return []
    
#%%
def getPTMs(row):
    
    modAA = re.findall('([A-Z])\[',row['Modified Peptide'])
    if not modAA:
        return 'NONE','NONE'
    modAAs = row['Modified Residues 3']


    dm= float(re.findall('\[(.*?)\]',row['Modified Peptide'])[0])
    
    
    
    #Check if the peptide N-terminus is modified
        #This returns the position of the first modified aa (with first aa position 1)
        #As indicated by the following [
        #It does not find the position of the second modification
    try:
        modpos = re.search('\[',row['Modified Peptide']).start()
    except:
        modpos = np.nan
    
    #Set tolerance to 25ppm 
    atol = row['Calibrated Observed M/Z'] *2.5*10**(-5)
    safebool = True
    if modAA:
        modAA = dict_123.get(modAA[0]) #Get 3 letter code of aa
    if modAAs:
        
        #Get filtered DFDM
        try:
            dfdmf = dfdm[modAAs]
        except:
            raise Exception(row)
        #Get index for SAAV which match the observed delta mass at the observed residue(s)    
        subindx = np.isclose(dfdmf,dm,atol=atol,rtol=0) 
        subindxloc = np.isclose(dfdm[modAA],dm,atol=atol,rtol=0) 
        #Get the index of common modifications which match the delta mass
        dangerindx = np.isclose(p.dfdangermods['Mass Shift'],dm,atol=atol,rtol=0)
        
        #Apply conditional filters for common modifications
        #Check if modified residues include those that match a common modification delta mass
        common_modified_aa = set(row['Modified Residues']).intersection(*p.dfdangermods.loc[dangerindx,'Modified Residues'])
        if common_modified_aa:
            safebool = False
        
        #If the modified residues don't explicitly apply to possible danger mods, check ambiguity
        if p.dfdangermods.loc[dangerindx,'Modified Residues'].str.contains('B').any():
           safebool = False #B is used for mods which apply to any residue
        #If the mod applies to a position
        if p.dfdangermods.loc[dangerindx,'Modified Position'].any():
            for possiblepos in p.dfdangermods.loc[dangerindx,'Modified Position']:
                if possiblepos == float(modpos):
                    safebool = False
            
        #Clear dangermod index if it failed to meet aa identity/position criteria
        if safebool: dangerindx = len(dangerindx)*[False]
           
        #Get all possible substitutions matching delta mass
        substitutions = ','.join([dfdmf.columns[col]+'->'+ dfdmf.index[row]
                                  for row, col in np.argwhere(subindx)])
        
        ptms = 'NONE'
        if not substitutions:
            substitutions = 'NONE'
        #Get assigned substitution per best localization
        if (dfdm.loc[subindxloc,modAA].any() & bool(p.dfdangermods[dangerindx]['Modification'].any())):
            ptms = modAA+'->'+','.join(dfdm.loc[subindxloc,modAA].index)+' or '+','.join(p.dfdangermods.loc[dangerindx,'Modification'].values)
        elif dfdm.loc[subindxloc,modAA].any():
            ptms = modAA+'->'+','.join(dfdm.loc[subindxloc,modAA].index)
        elif p.dfdangermods[dangerindx]['Modification'].any():
            ptms = ','.join(p.dfdangermods.loc[dangerindx,'Modification'].values)
        
    
        return ptms,substitutions
#%%
def FindTargets(dfin):#Given a peptide, mods and possible substitution sites, output possible sequences found
    df = dfin.copy()
    #Match origin and destination from ptms
    df['Origin'] = df['PTMs'].str.extract('(.{3})(?=->)')
    df['Destination'] = df['PTMs'].str.extract('(?<=->)(.{3})')
    maskIL = df['Destination'] == 'Leu'
    
    #Write substituted sequence
    df['Substituted Sequences'] = (df['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        df['Destination'].replace(p.dict_321) +
                        df['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0])
    #Write two sequences for I/L ambuguity
    dff = df[maskIL]
    dff.loc[:,'Substituted Sequences'] = (dff['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        'I' +
                        dff['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0] + ','+
                        dff['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        'L' +
                        dff['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0])
    
    df.loc[maskIL,'Substituted Sequences'] = dff['Substituted Sequences']
                        
    return df[['Origin','Destination','Substituted Sequences']]
#%%
def GetBaseBy(row):
    global bySeqFile,bySeq
    indx = (row['Sample'],row['Peptide'])
    try:
        result = bySeqFile.loc[indx]      #Try to get intensity matching file AND peptide
    except:
        try: 
            result = bySeq.loc[row['Peptide']]    #Try to get intensity matching only peptide
        except:
            #raise Exception('Failed to get groupby result for '+ str(indx))
            return None                    #No match returns none
    if result == 0:                        #If first attempt found 0 intensity, try just peptide match
        result = bySeq.loc[row['Peptide']]
    if result == 0:                        #If even matching only by peptide returns 0 intensity, return none
            return None
    return result
#%%
#Inputs
dfdm = pd.read_csv(r'C:\Users\taylo\anaconda3\Lib\tjl\params\dmMatrix.csv',index_col=0).astype(float)
if os.path.isdir(p.outputdir):
    print('Output directory already exists -- code will overwrite existing files with common names')
else: 
    os.mkdir(p.outputdir)
    print('Created output directory '+p.outputdir)


#Get Fragger filenames
filenames = pd.read_csv(os.path.join(p.inputdir,'filelist_ionquant.txt'),sep='\t')
#Drop spectra directory and flag column
filenames = filenames.loc[~filenames['flag'].str.contains('specdir'),'value']
#Drop psm.tsv
filenames = filenames.replace(regex=r'(\\.*)',value='')
filenames.drop_duplicates(inplace=True)
filenames.to_csv(os.path.join(p.outputdir,'filenames.csv'),index=False)

#Get raw filenames that correspond to Fragger filenames
fileSample = pd.read_csv(os.path.join(p.inputdir,"fragpipe-files.fp-manifest"),header=None ,sep='\t')
fileSample.columns=['Raw File','Sample','Replicate','Acquisition Type']
fileSample['Raw File'] = fileSample['Raw File'].apply(lambda x: os.path.split(x)[1])
#Remove .raw or .d extension
fileSample = fileSample.astype(str)
fileSample['Raw File'] = fileSample['Raw File'].str.replace('\.raw|\.d','')
fileSample = fileSample.set_index('Raw File')
#rawFileNames = '|'.join(fileSample['Raw File'])

columns = ['Spectrum','Peptide','Modified Peptide','Calibrated Observed Mass',
           'Intensity','Peptide Length','Assigned Modifications','Protein',
           'Hyperscore','PeptideProphet Probability','Retention',
           'Calibrated Observed M/Z','Number of Missed Cleavages']

#This information is useful for neighbor analysis, but was not present in earlier
#versions of file.psm output of FragPipe
columns += ['Prev AA','Protein Start','Protein End','MSFragger Localization']
#Comment this out if not using TIMSTOF data
columns =+ ['Ion Mobility']

#Pull out PSM info
df = pd.DataFrame()
for file in filenames:
    dfiter = pd.read_csv(os.path.join(p.inputdir,file,'psm.tsv'), sep='\t', usecols=columns)
    df = pd.concat([df,dfiter])
df.reset_index(drop=True,inplace=True)   

#Sometimes fragger hangs up on filenames with a .
#For runs where this needs to be renamed, I need to rename the spectrum data as well.
df['Spectrum'] = df['Spectrum'].str.replace(r'SALTY3.5_','SALTY3_5_')
fileSample.index = fileSample.index.str.replace(r'SALTY3.5_','SALTY3_5_')

#Replace Modified Peptide localization score [###] with delta mass 
df['Assigned Modifications'] = df['Assigned Modifications'].str.replace('\(57.021\d\)','')
df['Delta Mass'] = '['+df['Assigned Modifications'].str.extract('\((.*?)\)')+']'
df['Assigned Modifications'].fillna('NONE',inplace=True)
df['Modified Peptide'].fillna('NONE',inplace=True)
df['Delta Mass'].fillna('NONE',inplace=True)
df['MSFragger Localization'].fillna('NONE',inplace=True)

df['Modified Peptide'] = df.apply(lambda row: re.sub('\[(.*?)\]',row['Delta Mass'],
                                    row['Modified Peptide']), axis=1)


#Get possibly Modified Peptides
df['Modified Residues'] = df.apply(getModaa,axis=1)
df['Modified Residues 3'] = [[dict_123.get(aa) for aa in row if aa in dict_123] for row in df['Modified Residues']]
print('Localizations determined...')

#Find assigned modifications which match a sub
print('Naming modifications, this may take some time...')
now = datetime.datetime.now()
print(str(now))
#Get PTMs with best guess OR ambiguous MSFragger localization
df[['PTMs','All Possible Substitutions']] = df.apply(lambda x:getPTMs(x),
                                                axis=1,result_type='expand')

df.loc[:,'PTMs'] = df['PTMs'].fillna('NONE')


#Report binary Is Sub or Is Danger
print('Making binary filters...')
now = datetime.datetime.now()
print(str(now))

dfdecoy = df.loc[df['Protein'].str.contains('rev_')]
df = df.loc[~df['Protein'].str.contains('rev_')]
dfdecoy.to_csv(os.path.join(p.outputdir,'DiscoveryDecoyIons.csv'),index =False)

df['Is Sub'] = df['PTMs'].str.contains('->')
df['Is Danger'] = df['PTMs'].str.contains('or')

#Find peptides of any modification
df.loc[:,'Is Mod'] = ~(df['Modified Peptide'].str.contains('\['))

#Find Base peptides
modpeplist = df.loc[df['Is Sub'],'Peptide'].drop_duplicates()
df.loc[:,'Is Base'] = ~df['Is Sub']
df.loc[:,'Is Base'] = df['Is Base'] & df['Peptide'].isin(modpeplist)

#%%
#Get Modified Peptide

df[['Origin','Destination','Substituted Sequence']] = FindTargets(df)
df = p.ExpIL(df,'Substituted Sequence')

#%%
print('Getting Base Peptide Info...')
now = datetime.datetime.now()
print(str(now))

#Assign Sample names
df['Raw File'] = df['Spectrum'].str.extract('({})'.format('|'.join(fileSample.index)))
df['Sample'] = df['Raw File'].apply(lambda x: fileSample.loc[x,'Sample'])
df['Replicate'] = df['Raw File'].apply(lambda x: fileSample.loc[x,'Replicate'])

#Parse substituted peptides
dfsubs = df[df['Is Sub']]
dfsubs = dfsubs[~dfsubs['Is Danger']]
dfsubs.reset_index(inplace=True)

#Aggregate PSMs->Modified Sequences, by file, taking the max intensity
dfintensity = dfsubs.groupby(by=['Sample','Modified Peptide'])['Intensity'].max().to_frame()
dfintensity = dfintensity.merge(dfsubs.groupby(by=['Sample','Modified Peptide'])[['Peptide','Protein']].max(),
                                left_index=True,right_index=True)
dfintensity = dfintensity.reset_index()
#Same for tech rep level analysis
dftech = dfsubs.groupby(by=['Sample','Replicate','Modified Peptide'])['Intensity'].max().to_frame()
dftech = dftech.merge(dfsubs.groupby(by=['Sample','Replicate','Modified Peptide'])[['Peptide','Protein']].max(),
                                left_index=True,right_index=True)
dftech = dftech.reset_index()

#Aggregate PSMs->Modified Sequences for base peptides
dfbase = df[df['Is Base']]
bySeqFile = dfbase.groupby(by=['Sample','Peptide'])['Intensity'].max()
bySeq = dfbase.groupby(by=['Peptide'])['Intensity'].max()

#Get genomic cognate intensity
dfintensity['Base Intensity'] = dfintensity.apply(GetBaseBy,axis=1)

#Get Protein quantification
#Set common axis for dfintensity

dfintensity = dfintensity.set_index(['Protein','Sample'])
#Get protein quant
dfprot = pd.read_csv(os.path.join(p.inputdir,'combined_protein.tsv'),sep='\t')
protintcols = dfprot.columns[dfprot.columns.str.contains('Intensity')]
protintcols = ['Protein'] + protintcols[~protintcols.str.contains('Max')].to_list() 
dfprot = dfprot[protintcols].rename(columns=lambda x: x.replace(' Intensity', ''))
#Reshape and rename 
dfprot = dfprot.set_index('Protein').stack()
dfprot.name = 'Protein Intensity'
protbySample = dfprot.reset_index()
protbySample.columns = ['Protein','File','Protein Intensity']
protbySample['Sample'] = protbySample['File'].str.extract('(.*?)_')
protbySample = protbySample.groupby(by=['Protein','Sample'])['Protein Intensity'].max()
dfprot.index = dfprot.index.set_names(['Protein', 'Sample'])
dfintensity = dfintensity.merge(protbySample, how='left',left_index=True,right_index=True)
dftech = dftech.merge(dfintensity[['Base Intensity','Protein Intensity']],
                      how='left',left_on=['Protein','Sample'],right_index=True)
dfintensity = dfintensity.reset_index()
dfintensity = dfintensity[~dfintensity['Protein'].str.contains('(cont|rev_)')]
dftech = dftech[~dftech['Protein'].str.contains('(cont|rev_)')]
#Get normalized values

dfintensity['Fraction'] = dfintensity['Intensity']/(
    dfintensity['Base Intensity'] + dfintensity['Intensity'])
dfintensity['logFraction'] = np.log10(dfintensity['Fraction'].replace(0,np.nan))
dfintensity['Ratio'] = dfintensity['Intensity']/(dfintensity['Base Intensity'])
dfintensity['logRatio'] = np.log10(dfintensity['Ratio'].replace(0,np.nan))
dfintensity['Normalized'] = dfintensity['Ratio']/(dfintensity['Protein Intensity'])
dfintensity['logNormalized'] = np.log10(dfintensity['Normalized'].replace([0,np.inf,-np.inf],np.nan))

dftech['Fraction'] = dftech['Intensity']/(
    dftech['Base Intensity'] + dftech['Intensity'])
dftech['logFraction'] = np.log10(dftech['Fraction'].replace(0,np.nan))
dftech['Ratio'] = dftech['Intensity']/(dftech['Base Intensity'])
dftech['logRatio'] = np.log10(dftech['Ratio'].replace(0,np.nan))
dftech['Normalized'] = dftech['Ratio']/(dftech['Protein Intensity'])
dftech['logNormalized'] = np.log10(dftech['Normalized'].replace([0,np.inf,-np.inf],np.nan))

dfreproducability = dftech.groupby(by=['Sample','Modified Peptide'])[['logFraction','logRatio','logNormalized']].apply(
    lambda x: np.abs(np.std(x)/np.mean(x))*100)
#%%
#Create df for filtered only Ions, make general report
print('Writing tables...')
now = datetime.datetime.now()
print(str(now))

dftech.to_csv(os.path.join(p.outputdir,'SSP Quant Tech.csv'),index=False)
dfintensity.to_csv(os.path.join(p.outputdir,'SSP Quant.csv'),index=False)
dfsubs.to_csv(os.path.join(p.outputdir,'SSP PSM.csv'),index=False)
df.to_csv(os.path.join(p.outputdir,'AllPSMsAndFilters.csv'),index=False)


#%%
#Create global Substitution frequencies by sample
print('Making plots')
now = datetime.datetime.now()
print(str(now))

import seaborn as sns
from matplotlib import pyplot as plt

plt.figure()
sns.violinplot(data=dfintensity,x='Sample',y='logFraction')
plt.suptitle('Fraction normalized sub intensity')
plt.savefig(os.path.join(p.outputdir,'Fraction normalized sub intensity.png'))
plt.figure()
sns.violinplot(data=dfintensity,x='Sample',y='logRatio')
plt.suptitle('Ratio normalized sub intensity')
plt.savefig(os.path.join(p.outputdir,'Ratio normalized sub intensity.png'))
plt.figure()
sns.violinplot(data=dfintensity,x='Sample',y='logNormalized')
plt.suptitle('Protein normalized sub intensity')
plt.savefig(os.path.join(p.outputdir,'Protein normalized sub intensity.png'))

for col in dfreproducability.columns:
    plt.figure()
    plt.hist(dfreproducability[col],bins=50,range=(0,300))
    plt.xlabel(col+' CV%')
    plt.ylabel('# of substitutions')
    plt.savefig(os.path.join(p.outputdir,col+' CV.png'))

#%%
print('Done.')
now = datetime.datetime.now()
print(str(now))
import winsound
duration = 750  # milliseconds
freq = 500  # Hz
winsound.Beep(freq, duration)
winsound.Beep(freq, duration)
#Reset variables
#%reset -f

