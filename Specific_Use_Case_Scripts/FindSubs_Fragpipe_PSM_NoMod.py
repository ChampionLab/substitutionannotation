#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A derivative of FindSubs_Fragpipe_PSM.py, specifically for MSFragger searches that have the option "Report mass shift as a variable mod" set to "No"
FragPipe v22
Created on 2024.02.09
@author: taylorlundgren
"""

#%% Package Imports
import pandas as pd
import numpy as np
import os
import re 
import datetime
#from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from substitutionannotation.resources import assets as sA
from Bio import SeqIO
import warnings


#%%
def first_lowercase_index(string):
    if type(string) == str:
        lowercase = re.search("[a-z]", string)
        if lowercase:
         return int(lowercase.start())
        else:
         return 0
    else:
        return None

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
    
    if type(row['MSFragger Localization'])==str:
        
        unique_aa = [letter.upper() for letter in row['MSFragger Localization'] if letter.islower()]
        if not unique_aa:
            #For odd cases where Fragger has a mod but not a location
            #Return all aa in the sequence
            allaa = list(set(re.findall('[A-Z]',row['Modified Peptide'])))
            return  allaa
        else:
            return list(set(unique_aa))
    else:
        return []
    
#%%
        #<<< THIS IS THE CHANGE FOR THE NOMOD SETTING >>>
def getPTMs(row):
    if not type(row['MSFragger Localization'])==str:
        return None,None
    
    modAAs = row['Modified Residues 3']
    dm= row['Delta Mass']

    
    #Get position of first potentially modified aa, and its 3 letter code
    lowercase = re.search("[a-z]", row['MSFragger Localization'])
    if lowercase:
        modpos = lowercase.start()
        modAA = sA.dict_123ambiguous.get(row['MSFragger Localization'][modpos].upper())
    else:
        modAA = sA.dict_123ambiguous.get(row['MSFragger Localization'][0])
        modpos = 0
    
    
    #Set tolerance to 25ppm 
    atol = row['Calibrated Observed Mass'] *2.5*10**(-5) 
    isSub = True
    if modAAs:
        if 'U' in modAAs:
            warnings.warn(r'Modified amino acid U is not supported')
            return None,None

        #Get filtered DFDM
        try:
            dfdmf = dfdm[modAAs]
        except:
            raise Exception(row)
        #Get index for SAAV which match the observed delta mass at the observed residue(s)    
        subindx = np.isclose(dfdmf,dm,atol=atol,rtol=0) 
        subindxloc = np.isclose(dfdm[modAA],dm,atol=atol,rtol=0) 
        #Get the index of common modifications which match the delta mass
        dangerindx = np.isclose(sA.dfdangermods['Mass Shift'],dm,atol=atol,rtol=0)
        
        #Apply conditional filters for common modifications

        #Check if modified residues include those that match a common modification delta mass
        common_modified_aa = set(row['Modified Residues']).intersection(*sA.dfdangermods.loc[dangerindx,'Modified Residues'])
        if common_modified_aa:
            isSub = False
        
        #If the modified residues don't explicitly apply to possible danger mods, check ambiguity
        if sA.dfdangermods.loc[dangerindx,'Modified Residues'].str.contains('X').any():
           isSub = False #X is used for mods which apply to any residue
        #If the mod applies to a position
        if sA.dfdangermods.loc[dangerindx,'Modified Position'].any():
            for possiblepos in sA.dfdangermods.loc[dangerindx,'Modified Position']:
                if possiblepos == float(modpos):
                    isSub = False
            
        #Clear dangermod index if it failed to meet aa identity/position criteria
        if isSub: dangerindx = len(dangerindx)*[False]
           
        #Get all possible substitutions matching delta mass
        substitutions = ','.join([dfdmf.columns[col]+'->'+ dfdmf.index[row]
                                  for row, col in np.argwhere(subindx)])
        
        ptms = None
        if not substitutions:
            substitutions = None
        #Get assigned substitution per best localization
        if (dfdm.loc[subindxloc,modAA].any() & bool(sA.dfdangermods[dangerindx]['Modification'].any())):
            ptms = modAA+'->'+','.join(dfdm.loc[subindxloc,modAA].index)+' or '+','.join(sA.dfdangermods.loc[dangerindx,'Modification'].values)
        elif dfdm.loc[subindxloc,modAA].any():
            ptms = modAA+'->'+','.join(dfdm.loc[subindxloc,modAA].index)
        elif sA.dfdangermods[dangerindx]['Modification'].any():
            ptms = 'or '+','.join(sA.dfdangermods.loc[dangerindx,'Modification'].values)
        
    
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
                        df['Destination'].replace(sA.dict_321) +
                        df['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0])
    #Write two sequences for I/L ambuguity
    dff = df[maskIL]
    dff.loc[:,'Substituted Sequences'] = (dff['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        'I' +
                        dff['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0] + ','+
                        dff['Modified Peptide'].str.extract(r'(.*)(?=.\[)')[0] +
                        'L' +
                        dff['Modified Peptide'].str.extract(r'(?<=.\])(.*)')[0]).copy()
    
    df.loc[maskIL,'Substituted Sequences'] = dff['Substituted Sequences']
                        
    return df[['Origin','Destination','Substituted Sequences']]

#%%
def ExpIL(dfin,sequencecolumn): 
    """"Split list of possible sequences in provided
 df[sequencecolumn using , as delimiter, expanding the rows for each sequence."""
    df = dfin.copy()
    df[sequencecolumn] = df.loc[:,sequencecolumn].str.split(',')
    df = df.explode(sequencecolumn)
    return df
#%%
def assign_modified_sequence(row):
    if np.isnan(row['First Offset Index']):
        return None
    
    seq = row['Peptide']
    pos = int(row['First Offset Index']) +1
    out = seq[:pos] + '[' +  str(row['Mass Offset']) +']' + seq[pos:]
    return out


#%%

def find_codon(aa_pos,UniProtID):
    """
    Returns the genomic codon (3 letter string),
    for a given amino acid position (aa_pos, with 1 being the first aa),
    in a given protein (UniProtID)
    """

    sequence = dict_genome_sequences.get(UniProtID)
    if not sequence:
        return 'ERR_REF'
    #Note that I subtract 3 to align 0/1 indexing
    codon = dict_genome_sequences.get(UniProtID)[aa_pos*3-3:aa_pos*3]
    if not codon:
        return 'ERR_SEQ'
    return str(codon)

def find_codon_gene(aa_pos,gene):
    """
    Returns the genomic codon (3 letter string),
    for a given amino acid position (aa_pos, with 1 being the first aa),
    in a given protein (UniProtID)
    """

    sequence = dict_genome_sequences.get(gene)
    if not sequence:
        return 'ERR_REF'
    #Note that I subtract 3 to align 0/1 indexing
    codon = dict_genome_sequences.get(gene)[aa_pos*3-3:aa_pos*3]
    if not codon:
        return 'ERR_SEQ'
    return str(codon)

def get_mismatch_position(codon,possible_codon):
    mismatch_positions = ''
    if codon[0] != possible_codon[0]:
        mismatch_positions += "1"
    if codon[1] != possible_codon[1]:
        mismatch_positions += "2"
    if codon[2] != possible_codon[2]:
        mismatch_positions += "3"

    return mismatch_positions

#%%

def lowest_bp_mismatch(codon,aa):
    """
    Returns a tuple including the lowest number of mRNA/tRNA nucleotide base pair mismatches
    that result from a substitution from a given codon (codon),
    to a given amino acid (aa);
    Also including the positions of the mismatch for each combination of mismatches
    that still have the lowest number of mispaired nucleotides.
    """
    # Get the list of DNA codons for the given amino acid
    codons_for_aa = sA.codon_table_inverse[aa]
    if aa == 'L':   #Add Ile codons, because it is isobaric with Leu
        codons_for_aa += sA.codon_table_inverse['I']

    if not codons_for_aa:
        raise Exception(f'No codons found for {aa}')
    
    # Iterate through the codons for the amino acid
    mismatches = []
    mismatch_positions = []
    for possible_codon in codons_for_aa:
        # Calculate the character difference between the input codon and the codon from the table
        mismatches += [sum(c1 != c2 for c1, c2 in zip(codon, possible_codon))]
        try:
            mismatch_positions += [get_mismatch_position(codon,possible_codon)]
        except IndexError:
            raise Warning(f'Problems getting position of {codon} against {possible_codon}')
    #Get lowest number of mismatches. Sensibility check
    mismatch_min = np.min(mismatches)
    if mismatch_min == 0:
        warnings.warn('Cognate labeled as substitution?')
    if mismatch_min >3:
        raise Exception('Something went wrong counting bp mismatches...')

    #If multiple codons tie for lowest, get all possible mismatch positions
    best_mismatch_positions = []
    for imin in np.where(mismatches == mismatch_min)[0]:
        best_mismatch_positions += [mismatch_positions[imin]]

    return mismatch_min, best_mismatch_positions


def complement_dna(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complemented_sequence = ''.join(complement[base] for base in dna_sequence)
    return complemented_sequence

#%%
sA.timestamp('Importing Data...')

#Inputs
dfdm = sA.dfdm.copy()
if os.path.isdir(sA.outputdir):
    print('Output directory already exists -- code will overwrite existing files with common names')
else: 
    os.mkdir(sA.outputdir)
    print('Created output directory '+sA.outputdir)


#Get Fragger filenames
filenames = pd.read_csv(os.path.join(sA.inputdir,'filelist_ionquant.txt'),sep='\t')
#Drop spectra directory and flag column
filenames = filenames.loc[~filenames['flag'].str.contains('specdir'),'value']
#Drop psm.tsv
filenames = filenames.replace(regex=r'(\\.*)',value='')
filenames.drop_duplicates(inplace=True)
filenames.to_csv(os.path.join(sA.outputdir,'filenames.csv'),index=False)

#Get raw filenames that correspond to Fragger filenames
fileSample = pd.read_csv(os.path.join(sA.inputdir,"fragpipe-files.fp-manifest"),header=None ,sep='\t')
fileSample.columns=['Raw File','Sample','Replicate','Acquisition Type']
fileSample['Raw File'] = fileSample['Raw File'].apply(lambda x: os.path.basename(x))
#Remove .raw or .d extension
fileSample = fileSample.astype(str)
fileSample['Raw File'] = fileSample['Raw File'].str.replace('\.raw|\.d','',regex=True)
fileSample = fileSample.set_index('Raw File')

#Pull out PSM info
df = pd.DataFrame()
for file in filenames:
    dfiter = pd.read_csv(os.path.join(sA.inputdir,file,'psm.tsv'), sep='\t', low_memory=False)
    df = pd.concat([df,dfiter])
df.reset_index(drop=True,inplace=True)   


columns = ['Spectrum','Peptide','Modified Peptide','Delta Mass','Calibrated Observed Mass',
           'Intensity','Peptide Length','Assigned Modifications','Protein', 'Mapped Proteins',
           'Hyperscore','PeptideProphet Probability','Retention', 'Gene',
           'Calibrated Observed M/Z','Number of Missed Cleavages']

#This information is useful for neighbor analysis, but was not present in earlier
#versions of file.psm output of FragPipe
if (df.columns.str.contains('Prev AA').any() &
    df.columns.str.contains('Protein Start').any() &
    df.columns.str.contains('Protein End').any() &
    df.columns.str.contains('MSFragger Localization').any()):
    columns += ['Prev AA','Protein Start','Protein End','MSFragger Localization']
#Check if IM data is available
if df.columns.str.contains('Ion Mobility').any():
    columns += ['Ion Mobility']
#Remove unused columns
df = df[columns]

#Sometimes fragger hangs up on filenames with a .
#For runs where this needs to be renamed, I need to rename the spectrum data as well, otherwise there is discrepancy.
#Here is an example implementation for the data I used
df['Spectrum'] = df['Spectrum'].str.replace(r'SALTY3.5_','SALTY3_5_',regex=True)
fileSample.index = fileSample.index.str.replace(r'SALTY3.5_','SALTY3_5_',regex=True)

#Replace Modified Peptide localization score [###] with delta mass 
df['Assigned Modifications'] = df['Assigned Modifications'].str.replace('C\(57.021\d\)','',regex=True) #Removes fixed carbamidomethyl annotations
#df['Mass Offset'] = df['Assigned Modifications'].str.extract('\((.*?)\)').astype(float)
"""
#Resolve the odd case where a PTM is assigned a mass offset of 57.021 at a cysteine
#Note the mass offset is applied in addition to the delta mass of the fixed carbamidomethyl modification, these are not substitutions
doubleCysMask = (~df['Assigned Modifications'].isna() & df['Mass Offset'].isna())
print(f'Number of C+57.021 mods (after fixed carbamidomethylation): {doubleCysMask.sum()}')
df.loc[doubleCysMask,'Modified Peptide'] = df.loc[doubleCysMask,'Modified Peptide'].str.replace('\[.*\]','',regex=True)


df['Assigned Modifications'].fillna('NONE',inplace=True)
df['Modified Peptide'].fillna('NONE',inplace=True)
df['Mass Offset'] = df['Mass Offset'].astype(object)
#df['Mass Offset'].fillna('NONE',inplace=True)
df['MSFragger Localization'].fillna('NONE',inplace=True)
"""


#Get possibly Modified Peptides
df['Modified Residues'] = df.apply(getModaa,axis=1)
df['Modified Residues 3'] = [[sA.dict_123ambiguous.get(aa) for aa in row if aa in sA.dict_123ambiguous] for row in df['Modified Residues']]

df['Mass Offset'] = df['Delta Mass']
df['First Offset Index'] = df['MSFragger Localization'].apply(first_lowercase_index)
df['Modified Peptide'] = df.apply(assign_modified_sequence,axis=1)


#Find assigned modifications which match a sub
sA.timestamp('Naming modifications, this may take some time...')



# Prep the  progress bar
progress_bar = tqdm(total=len(df.index))
df[['PTMs','All Possible Substitutions']] = df.apply(lambda x:(getPTMs(x),progress_bar.update(1))[0],
                                                axis=1,result_type='expand')
# Close the progress bar
progress_bar.close()




#Report binary Is Sub or Is Danger
sA.timestamp('Making binary filters...')


dfdecoy = df.loc[df['Protein'].str.contains('rev_')]
df = df.loc[~df['Protein'].str.contains('rev_')]
dfdecoy.to_csv(os.path.join(sA.outputdir,'DiscoveryDecoyIons.csv'),index =False)

df['Is Sub'] = df['PTMs'].str.contains('->')
df['Is Sub'].fillna(False,inplace=True)
df['Is Danger'] = df['PTMs'].str.contains('or')
df['Is Danger'].fillna(False,inplace=True)

#Find peptides of any modification
df.loc[:,'Is Mod'] = df['Modified Peptide'].str.contains('\[')
df['Is Mod'].fillna(False,inplace=True)

#Find Base peptides
modpeplist = df.loc[df['Is Sub'],'Peptide'].drop_duplicates()
df.loc[:,'Is Base'] = ~df['Is Mod']
df.loc[:,'Is Base'] = df['Is Base'] & df['Peptide'].isin(modpeplist)

#%%
#Get Modified Peptide
sA.timestamp('Writing Substituted Sequences, then parsing...')
df[['Origin','Destination','Substituted Sequence']] = FindTargets(df)
df = ExpIL(df,'Substituted Sequence')

#%%

#Map Sample/replicate annotations to PSM data
df['Raw File'] = df['Spectrum'].str.extract('({})'.format('|'.join(fileSample.index)))
df['Sample'] = df['Raw File'].apply(lambda x: fileSample.loc[x,'Sample'])
df['Replicate'] = df['Raw File'].apply(lambda x: fileSample.loc[x,'Replicate'])

#Parse substituted peptides
dfsubs = df[df['Is Sub']]
dfsubs = dfsubs[~dfsubs['Is Danger']]                           #Remove ambiguous PTMs
dfsubs = dfsubs[~dfsubs['Protein'].str.contains('cont|rev_')]   #Remove contaminants and decoys
dfsubs['Intensity'] = dfsubs['Intensity'].replace(0,np.nan)     #Remove 0 intensity values for accurate mean/medians
dfsubs.reset_index(inplace=True)

#%%
sA.timestamp('Getting DNA codon information...')

#Get substituted position relative to the protein
dfsubs['Substitution Position'] =  dfsubs['Modified Peptide'].apply(
    lambda x:re.search('.\[',x).start()) + dfsubs['Protein Start']
dfsubs['UniProt ID'] = dfsubs['Protein'].str.extract(r'\|(.*)\|')

# Load the genome sequence from the GenBank file
genome_record = SeqIO.read(sA.path_to_GenBank, "genbank")
genome_sequence = genome_record.seq

#Get a dictionary with UniProt ID as the key and genomic sequence as the value
dict_genome_sequences = {}

     
#db reference based alignment
for feature in genome_record.features:
    try:
        for ref in feature.qualifiers['db_xref']:
            if ref.find('UniProt') == 0:
                uniprot = re.findall('UniProt.*:(.*)',ref)[0]
                start = feature.location.start
                end = feature.location.end
                if dict_genome_sequences.get(uniprot):
                    warnings.warn(f'Duplicate {uniprot} entires')
                else:   # genome_sequence is the top strand, we need to orient proteins translated in order of bottom strand dna
                    if feature.strand == 1:         #Top strand
                        dict_genome_sequences[uniprot] = genome_sequence[start:end]
                    if feature.strand == -1:   #Bottom Strand, -1 offsets the [,) inclusivity of python slicing
                        dict_genome_sequences[uniprot] = complement_dna(genome_sequence[end-1:start-1:-1])
    
    except KeyError:
        continue    #Not all the genome features will have db_xref
dfsubs['Codon'] = dfsubs.apply(lambda row: find_codon(row['Substitution Position'],row['UniProt ID']),axis=1)

"""
#gene name based alignment
#Note this is still problematic - many duplicate gene features by name
for feature in genome_record.features:
    try:
        gene = feature.qualifiers['gene'][0]
        start = feature.location.start
        end = feature.location.end
        if dict_genome_sequences.get(gene):
            warnings.warn(f'Duplicate {gene} entires')
        else:   # genome_sequence is the top strand, we need to orient proteins translated in order of bottom strand dna
            if feature.strand == 1:         #Top strand
                dict_genome_sequences[gene] = genome_sequence[start:end]
            if feature.strand == -1:   #Bottom Strand, -1 offsets the [,) inclusivity of python slicing
                dict_genome_sequences[gene] = complement_dna(genome_sequence[end-1:start-1:-1])
    
    except KeyError:
        continue    #Not all the genome features will have gene

dfsubs['Codon'] = dfsubs.apply(lambda row: find_codon_gene(row['Substitution Position'],row['Gene']),axis=1)
"""

dfsubs[['# BP Mismatches','BP Mismatch Positions']] = dfsubs.apply(lambda row: lowest_bp_mismatch(row['Codon'],sA.dict_321.get(row['Destination'])),axis=1,result_type='expand')

#%%

sA.timestamp('Getting quantifications...')

#Prep base peptides
dfbase = df[df['Is Base']]

#Get Protein quantification
    #Get protein quant from MaxLFQ Intensity
dfprot = pd.read_csv(os.path.join(sA.inputdir,'combined_protein.tsv'),sep='\t')
protintcols = dfprot.columns[dfprot.columns.str.contains('Intensity')]  
protintcols = protintcols[protintcols.str.contains('MaxLFQ',regex=True)].to_list() + ['Protein']
dfprot = dfprot[protintcols].rename(columns=lambda x: x.replace(' MaxLFQ Intensity', ''))
    #Reshape and rename, annotate sample/replicate information 
dfprot = dfprot.set_index('Protein').stack()
dfprot.name = 'Protein Intensity'
protbySample = dfprot.reset_index()
protbySample.columns = ['Protein','File','Protein Intensity']

protbySample['Sample'] = protbySample['File'].str.extract('(.*?)_')
protbySample['Replicate'] = protbySample['File'].str.extract('_(\d*)')

#Handle mismatch of sample annotations 
# Including from older versions of fragger without fragpipe manifest,
# or more generally when 'sample' was set to file name but here I redefine the samples

#For single biological replicates / no replicates input into FragPipe
if protbySample['File'].str.contains('_').any() == False:
    protbySample['Sample'] = protbySample['File']
    protbySample['Replicate'] = 1
    protbySample['Replicate'] = protbySample['Replicate'].astype(object)

# This probably only applies to my old stuff
if not protbySample['Sample'].isin(fileSample['Sample']).any():
    
    try:    #This is what I used for old data 
        protbySample['Sample'] = protbySample.apply(lambda row: fileSample.loc[fileSample.index==row['File'],'Sample'].values.item(),axis=1)
        protbySample['Replicate'] = protbySample.apply(lambda row: fileSample.loc[fileSample.index==row['File'],'Replicate'].values.item(),axis=1)
    except ValueError: #If this logic doesn't align sample annotation
        try:    #This is for a specific case where an _ in the sample name messed things up
          protbySample['Sample'] = protbySample['File'].str.split('_').apply(lambda x: '_'.join(x[:-1]))
          protbySample['Replicate'] = protbySample['File'].str.split('_').apply(lambda x: x[-1])
            #Check that all extracted values align
          if not protbySample['Sample'].isin(fileSample['Sample']).all():
            raise Exception('Could not align sample annotation with combined_protein.csv')
        except:
            raise Exception('Could not align sample annotation with combined_protein.csv')

protbySample['Protein Intensity'] = protbySample['Protein Intensity'].replace(0,np.nan) 
protbySample.drop(columns='File',inplace=True)


#Start by aggregating technical replicates ('File') into summed biological replicate data
byBio = dfsubs.groupby(by=['Sample','Replicate','Modified Peptide'])

    #Protein and cognate Peptide columns all share the same value, so max is just an easy tool
bioQuant = byBio['Protein'].max().to_frame()
bioQuant['Peptide'] = byBio['Peptide'].max()
bioQuant['Gene'] = byBio['Gene'].max()

    #Sum intensity of technical replicates within each sample, biological replicate for each modified peptide
bioQuant['Intensity'] = byBio['Intensity'].sum().replace(0,np.nan)
    #Do the same for genomic cognate peptides
baseBioQuant = dfbase.groupby(by=['Sample','Replicate','Peptide'])['Intensity'].sum().replace(0,np.nan)
baseBioQuant.rename('Base Intensity',inplace=True)
    #Count the number of non-0 intensity technical replicates
    #Note I need to aggregate by File as well, so I don't count multiple PSMs from each individual run/'File'
byFile = dfsubs.groupby(by=['Sample','Replicate','Raw File','Modified Peptide'])['Intensity'].max().dropna()
bioQuant['Number of Technical Replicates'] = byFile.reset_index().groupby(by=['Sample','Replicate','Modified Peptide'])['Intensity'].apply(len)
    #Merge in genomic cognate and protein information last (indexing reasons)
bioQuant = bioQuant.reset_index().merge(baseBioQuant,on=['Sample','Replicate','Peptide'])
bioQuant = bioQuant.merge(protbySample,how='left',on=['Protein','Sample','Replicate'])


#Aggregate biological replicate data to sample data
bySample = bioQuant.groupby(by=['Sample','Modified Peptide'])
    #Average replicate intensities for biological replicate intensity
sampleQuant = bySample['Intensity'].mean().replace(0,np.nan).to_frame()
sampleQuant['Base Intensity'] = bySample['Base Intensity'].mean().replace(0,np.nan)
sampleQuant['Protein Intensity'] = bySample['Protein Intensity'].mean().replace(0,np.nan)

    #Protein and cognate Peptide columns all share the same value, so max is just an easy tool
sampleQuant['Protein'] = bySample['Protein'].max()
sampleQuant['Peptide'] = bySample['Peptide'].max()
sampleQuant['Peptide'] = bySample['Gene'].max()

    #Count the number of non-0 intensity biological replicates
sampleQuant['Number of Biological Replicates'] = bioQuant.reset_index().groupby(by=['Sample','Modified Peptide'])['Intensity'].apply(lambda x: len(x.replace(0,np.nan).dropna()))

#Calculate normalized values
sampleQuant['Ratio'] = sampleQuant['Intensity']/sampleQuant['Base Intensity']
sampleQuant['Fraction'] = sampleQuant['Intensity']/(
    sampleQuant['Base Intensity'] + sampleQuant['Intensity'])
sampleQuant['Protein Normalized'] = sampleQuant['Intensity']/sampleQuant['Protein Intensity']

bioQuant['Ratio'] = bioQuant['Intensity']/bioQuant['Base Intensity']
bioQuant['Fraction'] = bioQuant['Intensity']/(
    bioQuant['Base Intensity'] + bioQuant['Intensity'])
bioQuant['Protein Normalized'] = bioQuant['Intensity']/bioQuant['Protein Intensity']


#Calculate log10 Values
sampleQuant['logRatio'] = np.log10(sampleQuant['Ratio'])
sampleQuant['logFraction'] = np.log10(sampleQuant['Fraction'])
sampleQuant['logProtein Normalized'] = np.log10(sampleQuant['Protein Normalized'])

bioQuant['logRatio'] = np.log10(bioQuant['Ratio'])
bioQuant['logFraction'] = np.log10(bioQuant['Fraction'])
bioQuant['logProtein Normalized'] = np.log10(bioQuant['Protein Normalized'])

#Calculate statistical metrics
sampleQuant['Sub SD'] = bySample['Intensity'].std()
sampleQuant['Cognate SD'] = bySample['Base Intensity'].std()
sampleQuant['Protein CV'] = bySample['Protein Intensity'].std()
sampleQuant['Sub CV'] = sampleQuant['Sub SD']/sampleQuant['Intensity']
sampleQuant['Cognate CV'] = sampleQuant['Cognate SD']/sampleQuant['Base Intensity']
sampleQuant['Ratio CV'] = np.sqrt(sampleQuant['Sub CV']**2 + sampleQuant['Cognate CV']**2)
#sampleQuant['Fraction CV'] =  Maybe I'll add this later
sampleQuant['Protein Normalized CV'] = np.sqrt(sampleQuant['Sub CV']**2 + sampleQuant['Protein CV']**2)

sampleQuant.reset_index(inplace=True)

try: protbySample['Replicate'] = protbySample['Replicate'].astype(int)
except: raise Warning('Failed to interpret replicate annotation')

if protbySample['Replicate'].max() > 1:
    quantiles = pd.DataFrame()
    quantiles['Sub CV'] = np.quantile(sampleQuant['Sub CV'].dropna(),.95)
    quantiles['Cognate CV'] = np.quantile(sampleQuant['Cognate CV'].dropna(),.95)
    quantiles['Ratio CV'] = np.quantile(sampleQuant['Ratio CV'].dropna(),.95)
    quantiles['Protein CV'] = np.quantile(sampleQuant['Protein CV'].dropna(),.95)
    quantiles.to_csv(os.path.join(sA.outputdir,'CV95 Quantiles.csv'),index=False)


#%%
#Write tables
print('Writing tables...')
now = datetime.datetime.now()
print(str(now))

bioQuant.to_csv(os.path.join(sA.outputdir,'SSP Quant Replicates.csv'),index=False)
sampleQuant.to_csv(os.path.join(sA.outputdir,'SSP Quant.csv'),index=False)
dfsubs.to_csv(os.path.join(sA.outputdir,'SSP PSM.csv'),index=False)
baseBioQuant.reset_index().to_csv(os.path.join(sA.outputdir,'Base Peptides.csv'))
df.to_csv(os.path.join(sA.outputdir,'AllPSMsAndFilters.csv'),index=False)



#%%
#Plot global Substitution frequencies by sample
print('Making plots')
now = datetime.datetime.now()
print(str(now))

import seaborn as sns
from matplotlib import pyplot as plt

plt.figure()
sns.violinplot(data=sampleQuant,x='Sample',y='logFraction',cut=0)
plt.suptitle('Fraction normalized sub intensity')
plt.savefig(os.path.join(sA.outputdir,'Fraction normalized sub intensity.png'))
plt.figure()
sns.violinplot(data=sampleQuant,x='Sample',y='logRatio',cut=0)
plt.suptitle('Ratio normalized sub intensity')
plt.savefig(os.path.join(sA.outputdir,'Ratio normalized sub intensity.png'))
plt.figure()
sns.violinplot(data=sampleQuant,x='Sample',y='logProtein Normalized',cut=0)
plt.suptitle('Protein normalized sub ratio')
plt.savefig(os.path.join(sA.outputdir,'Protein normalized sub ratio.png'))



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


# %%
