#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A script to annotate the output of Fragpipe mass-offset search, identifying peptides with an amino acid substitution.
MSFragger must have the delta mass of each substitution and other potentially interfering post-translational modifications included in the delta mass option.
MSFragger must have the option to report mass-offset as a variable modification was set to 1 (“Yes - and remove delta mass”). 
    This impacts how the mass accuracy is referenced in determining if a possible PTM is mass-ambiguous with a PSM.
    This also affects how IonQuant will quantify the peptides.
This script only filters based on PTM type and isobaric PTMs, and does not include other common filters such as the number of replicate observations.
If you adjust the common PTM table, you MUST specify the delta mass AND the modified residues in the table. Otherwise, they will not be assigned and 
prevent erroneous substitution assignments. You may use 'X' as a wildcard for any residue. Specifying a modified PEPTIDE position is optional.
This script also assumes that decoy hits are annotated as _rev (Fragger default) or omitted.


Created on Thu Feb 25 11:43:35 2021
@author: taylorlundgren
"""

#%% Package Imports
import pandas as pd
import numpy as np
import os
import re 
import datetime
#from concurrent.futures import ProcessPoolExecutor, as_completed           Multithreaded moved to FindSubs_Fragpipe_multi.py
from tqdm import tqdm
from substitutionannotation.resources import assets as sA
from Bio import SeqIO
import warnings
import ast

#%%
complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
dict_to_mrna = {'T':'U','A':'A','C':'C','G':'G'}


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
        
        unique_aa = [letter for letter in row['MSFragger Localization'] if letter.islower()]
        unique_aa = [letter.upper() for letter in unique_aa]
        if not unique_aa:
            #For odd cases where Fragger has a mod but not a location
            #Return all aa in the sequence
            allaa = list(set(re.findall('[A-Z]',row['Modified Peptide'])))
            return  allaa
        else:
            return unique_aa
    else:
        return []
    
#%%
def getPTMs(row):
    
    modAA = re.findall('([A-Z])\[',row['Modified Peptide'])
    if not modAA:
        return 'NONE','NONE'
    modAAs = row['Modified Residues 3']

    #With current settings for IonQuant to work, MSFragger adds the delta mass to the
    #calculated peptide mass. The true delta mass from the base peptide is therefore
    #the sum of the reported delta mass and the reported mass offset
    #Note, the reported mass offset is NOT always the offset that minimizes the 
    #reported delta mass... why? Probably a MS2 thing but I haven't looked into it
    dm= row['Delta Mass'] + row['Mass Offset']

    
    #Check if the peptide N-terminus is modified
        #This returns the position of the first modified aa (with first aa position 1)
        #As indicated by the following [
        #It does not find the position of a second modification
    try:
        modpos = re.search(r'\[',row['Modified Peptide']).start()
    except:
        modpos = np.nan
    try:        #Change the modpos to -1 if the last aa is modified.
                #Note this requires you remove any fixed/variable modified annotations [] from the Modified Peptide column.
        if re.search(r'\]',row['Modified Peptide']).end() == len(row['Modified Peptide']):
            modpos = -1
    except:
        modpos = np.nan
    
    #Set tolerance to 25ppm 
    atol = row['Calibrated Observed Mass'] *2.5*10**(-5) 
    isSub = True
    if modAAs:
        if 'U' in modAAs:
            warnings.warn(r'Modified amino acid U is not supported')
            return 'UNK','UNK'
        modAA = sA.dict_123ambiguous.get(modAA[0]) #Get 3 letter code of aa

        #Get filtered DFDM
        try:
            dfdmf = dfdm[modAAs]
        except:
            raise Exception(f'Could not align dfdm with {row}')
        #Get index for SAAV which match the observed delta mass at the observed residue(s)    
        subindx = np.isclose(dfdmf,dm,atol=atol,rtol=0) 
        subindxloc = np.isclose(dfdm[modAA],dm,atol=atol,rtol=0) 
        #Get the index of common modifications which match the delta mass
        dangerindx = np.isclose(sA.dfdangermods['Mass Shift'],dm,atol=atol,rtol=0)
        
        #Apply conditional filters for common modifications

        #Check if modified residues include those that match a common modification delta mass
        danger_residues = set()     #Get residues that are modified by common mods
        for item in sA.dfdangermods.loc[dangerindx,'Modified Residues']:
            mod_residues = set(item.split(','))
            danger_residues.update(mod_residues)
            #Find common residues between the modified residues and the common modifications
        common_modified_aa = set(row['Modified Residues']).intersection(danger_residues)
        if common_modified_aa:
            isSub = False
        
        #If the modified residues don't explicitly apply to possible danger mods, check ambiguity
        if 'X' in danger_residues:
           isSub = False #X is used for mods which apply to any residue
        
        #If the mod applies to a position
        if sA.dfdangermods.loc[dangerindx,'Modified Position'].any():
            for possiblepos in sA.dfdangermods.loc[dangerindx,'Modified Position']:
                if possiblepos == float(modpos):
                    isSub = False
            
        #Clear dangermod index if it failed to meet aa identity/position criteria
        if isSub: 
            dangerindx = len(dangerindx)*[False]
           
        #Get all possible substitutions matching delta mass
        substitutions = ','.join([dfdmf.columns[col]+'->'+ dfdmf.index[row]
                                  for row, col in np.argwhere(subindx)])
        
        ptms = 'UNK'
        if not substitutions:
            substitutions = 'NONE'
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
    df['Destination'] = df['PTMs'].str.extract(r'(?<=->)(.{3})')
    maskIL = df['Destination'] == 'Xle'
    
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
                        
    return df[['Destination','Substituted Sequences']]

#%%
def ExpIL(dfin,sequencecolumn): 
    """"Split list of possible sequences in provided
 df[sequencecolumn using , as delimiter, expanding the rows for each sequence."""
    df = dfin.copy()
    df[sequencecolumn] = df.loc[:,sequencecolumn].str.split(',')
    df = df.explode(sequencecolumn)
    return df


#%%
def get_mismatch_position(codon,possible_codon):
    mismatch_positions = ''
    if codon[0] != possible_codon[0]:
        mismatch_positions += "1"
    if codon[1] != possible_codon[1]:
        mismatch_positions += "2"
    if codon[2] != possible_codon[2]:
        mismatch_positions += "3"

    return mismatch_positions

def lowest_bp_mismatch(row):
    """
    Returns a tuple including the lowest number of mRNA/tRNA nucleotide base pair mismatches
    that result from a substitution from a given codon (codon),
    to a given amino acid (aa);
    Also including the positions of the mismatch for each combination of mismatches
    that still have the lowest number of mispaired nucleotides.
    """
    progress_bar.update()
    codon = row['Codon']
    aa = sA.dict_321.get(row['Destination'])
    row_idx = row.name
    #Some codons were not properly identified (for different reasons). Return none for these codons...
    if codon == 'Err':
        return None,None,None
    if codon == 'ERR':
        return None,None,None

    # Get the list of DNA codons for the given amino acid
    codons_for_aa = sA.codon_table_inverse[aa].copy()
    if aa == 'L':   #Add Ile codons, because it is isobaric with Leu
        codons_for_aa += sA.codon_table_inverse['I'].copy()

    if not codons_for_aa:
        raise Exception(f'No codons found for {aa} row {row_idx}')
    
    # Iterate through the codons for the amino acid
    mismatches = []
    nmismatches = []
    mismatch_positions = []
    for possible_codon in codons_for_aa:
        # Calculate the character difference between the input codon and the codon from the table
        nmismatches.append(sum([c1 != c2 for c1, c2 in zip(codon, possible_codon)]))
        mismatches.append(str.join('',[dict_to_mrna[c1] if c1 != c2 else '' for c1, c2 in zip(codon, possible_codon)]))
        try:
            mismatch_positions += [get_mismatch_position(codon,possible_codon)]
        except IndexError:
            raise Warning(f'Problems getting position of {codon} against {possible_codon} row {row_idx}')
    #Get lowest number of mismatches. Sensibility check
    mismatch_min = np.min(nmismatches)
    if mismatch_min == 0:
        if row['Substitution Position'] == 1:
            warnings.warn('Detected non-Met aa at protein N-terminus!')
        else:
            warnings.warn(f'Cognate labeled as substitution? Row {row_idx}')
    if mismatch_min >3:
        raise Exception(f'Something went wrong counting bp mismatches... row {row_idx}')

    #If multiple codons tie for lowest, get all possible mismatch positions
    best_mismatch_positions = []
    best_mismatches = []
    for imin in np.where(nmismatches == mismatch_min)[0]:
        best_mismatch_positions += [mismatch_positions[imin]]
        best_mismatches += [mismatches[imin]]

    return mismatch_min, best_mismatch_positions, best_mismatches


def complement_dna(dna_sequence):
    complemented_sequence = ''.join(complement[base] for base in dna_sequence)
    return complemented_sequence

#%%
sA.timestamp('Importing Data...')

#Inputs
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
fileSample['Raw File'] = fileSample['Raw File'].str.replace(r'\.raw|\.d','',regex=True)
fileSample = fileSample.set_index('Raw File')

#Identify if carbamidomethylation was used as a variable mod. If not, assume it was used as a fixed mod.
has_varmod = False
active = False
pattern = re.compile(r'#?\s*variable_mod_\d{2}\s*=\s*57\.02146\s*C\s*1')
with open(os.path.join(sA.inputdir,'fragger.params'), 'r') as file:
    for line in file:
        if pattern.match(line):
            has_varmod = True
            if not line.strip().startswith('#'):
                active = True
            break  # Stop after the first match, remove this if you need to check all lines

if has_varmod & active:
    dfdm = sA.dfdm_varcys.copy()   #DFDM reflects genomic, unmodified Cys and substituted, camCys
else:
    dfdm = sA.dfdm.copy()           #DFDM reflects camCys both genomically and substituted

#Pull out PSM info
df = pd.DataFrame()
for file in filenames:
    dfiter = pd.read_csv(os.path.join(sA.inputdir,file,'psm.tsv'), sep='\t', low_memory=False)
    df = pd.concat([df,dfiter])
df.reset_index(drop=True,inplace=True)   


columns = ['Spectrum','Peptide','Modified Peptide','Delta Mass','Calibrated Observed Mass',
           'Intensity','Peptide Length','Assigned Modifications','Protein', 'Mapped Proteins','Protein ID',
           'Hyperscore','Retention', 'Gene',
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
#Fragpipe v22 changed some outputs for PSM scoring...
if df.columns.str.contains('PeptideProphet Probability').any():
    columns += ['PeptideProphet Probability']
elif df.columns.str.contains('SpectralSim').any():
    columns += ['SpectralSim','RTScore','Probability']
#Remove unused columns
df = df[columns]
#Parse gene information from Protein ID - this may be in Fragger's 'Gene' column if your input fasta is formatted right
df['Gene'] = df['Protein ID'].str.extract(r'gene=([^]]+)')

#Sometimes fragger hangs up on filenames with a .
#For runs where this needs to be renamed, I need to rename the spectrum data as well, otherwise there is discrepancy.
#Here is an example implementation for the data I used
df['Spectrum'] = df['Spectrum'].str.replace(r'SALTY3.5_','SALTY3_5_',regex=True)
fileSample.index = fileSample.index.str.replace(r'SALTY3.5_','SALTY3_5_',regex=True)


def insert_offsets(modified_peptide,assigned_modifications):
    # Ensure assigned_modifications exists
    if not isinstance(assigned_modifications, str):
        return modified_peptide  # Return original modified_peptide if assigned_modifications is not a string
        #Get offsets in alphabetical order (i.e. 12A[14] would be before 6A[14])
    
    #Some mass-offsets are terminal specific, and have...
        # Assigned Modifications annotation as 'N-term(mass offset number)' 
        # Modified Peptide with an additional n[mass offset number]PEPTIDESEQUENCE instead of P[mass offset number]EPTIDESEQUENCE
        #as of MSFragger v22.0. Let's adjust this to be compatible with the pipeline...
    assigned_modifications = re.sub(r'N-term','1X',assigned_modifications)    
    if modified_peptide[0] == 'n':
        modified_peptide = re.sub(r'(n\[\d+\])','',modified_peptide,1)
        modified_peptide = modified_peptide[0]+'[100]'+modified_peptide[1:] 
         
    offsets = re.findall(r'\((.*?)\)',assigned_modifications)
    positions = re.findall(r'(\d+)[A-Z]\(',assigned_modifications)
        #Tuple original index with position value, then sort by position value
    sorted_positions = [(i,int(x)) for i,x in enumerate(positions)]
    sorted_positions.sort(key=lambda x:x[1])
        #Use ordered original index to order offsets by position value
    sorted_offsets = [offsets[i] for i,_ in sorted_positions]
        #For each offset/modification, replace the annotation in 'Modified Peptide'
            #This makes 'Modified Peptide' a representation of each unique sequence, used a lot downstream
    for offset in sorted_offsets:
        offset = '('+offset+')'     #Change [] to () so we always overwrite one new value
        modified_peptide = re.sub(r'\[(\d+)\]',offset,modified_peptide,1)
    modified_peptide = re.sub(r'\(','[',modified_peptide)        #Change back to [] for downstream annotations
    modified_peptide = re.sub(r'\)',']',modified_peptide)
    return modified_peptide

df['Modified Peptide'] = df.apply(lambda row: insert_offsets(row['Modified Peptide'],row['Assigned Modifications']),axis=1)
df['Modified Peptide'] = df['Modified Peptide'].str.replace(r'C\[57.021\d\]','C',regex=True)                  #Removes carbamidomethyl annotations
df['Assigned Modifications'] = df['Assigned Modifications'].str.replace(r'C\(57.021\d\)','',regex=True)      #Removes carbamidomethyl annotations
df['Mass Offset'] = df['Assigned Modifications'].str.extract(r'\((.*?)\)').astype(float)
df['Assigned Modifications'] = df['Assigned Modifications'].fillna('NONE')
df['Modified Peptide'] = df['Modified Peptide'].fillna('NONE')
df['Mass Offset'] = df['Mass Offset'].astype(float)

#df['Mass Offset'].fillna('NONE',inplace=True)
df['MSFragger Localization'] = df['MSFragger Localization'].fillna('NONE')


#Get possibly Modified Peptides
df['Modified Residues'] = df.apply(getModaa,axis=1)
df['Modified Residues 3'] = [[sA.dict_123ambiguous.get(aa) for aa in row if aa in sA.dict_123ambiguous] for row in df['Modified Residues']]

#Find assigned modifications which match a sub
sA.timestamp('Naming modifications, this may take some time...')

# Prep the  progress bar
progress_bar = tqdm(total=len(df.index))
df[['PTMs','All Possible Substitutions']] = df.apply(lambda x:(getPTMs(x),progress_bar.update(1))[0],
                                                axis=1,result_type='expand')
# Close the progress bar
progress_bar.close()


df.loc[:,'PTMs'] = df['PTMs'].fillna('NONE')


#Report binary Is Sub or Is Danger
sA.timestamp('Making binary filters...')


dfdecoy = df.loc[df['Protein'].str.contains('rev_')]
df = df.loc[~df['Protein'].str.contains('rev_')]
df = df.loc[~df['Protein'].str.contains('cont_')]
dfdecoy.to_csv(os.path.join(sA.outputdir,'DiscoveryDecoyIons.csv'),index =False)

df['Is Sub'] = df['PTMs'].str.contains('->')
df['Is Danger'] = df['PTMs'].str.contains('or')
df.loc[df['Is Danger'], 'Is Sub'] = False           #Remove 'Is Sub' label from ambiguous substitutions

#Find peptides of any modification
df.loc[:,'Is Mod'] = df['Modified Peptide'].str.contains(r'\[')

#Find Base peptides
modpeplist = df.loc[df['Is Sub'],'Peptide'].drop_duplicates()
df.loc[:,'Is Base'] = ~df['Is Mod']
df.loc[:,'Is Base'] = df['Is Base'] & df['Peptide'].isin(modpeplist)

#%%
#Get Modified Peptide
sA.timestamp('Writing Substituted Sequences, then parsing...')
df[['Destination','Substituted Sequence']] = FindTargets(df)

def find_origin(modpep):
    if re.findall(r'([A-Z])\[',modpep):
        return  sA.dict_123.get(re.findall(r'([A-Z])\[',modpep)[0])
    else:
        return None

df['Origin'] = df['Modified Peptide'].apply(find_origin)

#Xle Substituted sequence is written as comma seperated with both options.
    #I chose to double the PSM, with each unique sequence assigned.
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

#This assumes that you use the same accessions for your protein and dna FASTA
genome_index = SeqIO.index(sA.path_dna_fasta,format='fasta')
#Get substituted position relative to the protein
dfsubs['Substitution Position'] =  dfsubs['Modified Peptide'].apply(
    lambda x:re.search(r'.\[',x).start()) + dfsubs['Protein Start']

prot_no_gene = []
bad_codon_index = []
def get_codon(prot,aa_position,origin,mapped_proteins,idx):
    """
    A function to find the codon corresponding to the location within a protein where a substitution was identified.
    Designed to be used row-wise in a dataframe with substitution PSMs on each row.
    prot = protein ID
    aa_position = position relative to the protein, input expects count to start at 1
    origin = the aa expected by the genomic sequence in the protein FASTA. Note this sometimes does not align with translation of the genomic FASTA,
        such as the first aa as Met in the protein sequence but being encoded by a Val or Leu codon.
    idx = the index of the row being operated on, for code troubleshooting purposes
    """
    #Get the codon
    gene = prot.replace('_prot_','_cds_')
    try:
        dna = genome_index.get(gene).seq
    except AttributeError:
        prot_no_gene.append(prot)
        return 'Err'
    start = (aa_position-1)*3
    codon = dna[start:start+3]
    #Logical checks on detected codon
    if str(codon.translate()) != sA.dict_321.get(origin):       #Does this match protein sequence in search database?
        if (origin == 'Met') & (aa_position == 1):                #My protein database coerces the first aa to Met, but sometimes first aa is not AUG
            return(str(codon))
        if (str(codon.translate()) in ['I','L']) & (sA.dict_321.get(origin) in ['I','L']): #Check for I/L ambiguity
            return(str(codon))
        
        if type(mapped_proteins) == str:                                     #Check ambiguous proteins in protein group for a better match
            mapped_proteins = mapped_proteins.split(',')                     
            for other_prot in mapped_proteins: 
                other_gene = other_prot.replace('_prot_','_cds_')                    
                try:
                    other_dna = genome_index.get(other_gene).seq
                except AttributeError:
                    continue
                other_codon = other_dna[start:start+3]
                if other_codon.translate() == sA.dict_321.get(origin):
                    print('Codon recovered from mapped proteins')
                    return str(other_codon)

        bad_codon_index.append(idx)
        return 'ERR'
    return str(codon)

dfsubs['Codon'] = dfsubs.apply(lambda row: get_codon(row['Protein'],row['Substitution Position'],row['Origin'],row['Mapped Proteins'],row.name),axis=1)
if len(prot_no_gene)>0:
    print(f'Failed to get DNA sequence for {np.unique(prot_no_gene)}')
    print('These proteins will be removed from reported tables. Please update the reference FASTA file to keep these observations.')
print(f'Got the wrong codon for {len(bad_codon_index)} PSMs, index stored in bad_codon_index list.')

#Filter proteins not found in genomic FASTA
df = df[~df['Protein'].isin(prot_no_gene)]
dfsubs = dfsubs[~dfsubs['Protein'].isin(prot_no_gene)]

print('Looking up codons in tables')
#Import lookup tables
table_mismatch_count = pd.read_csv(os.path.join(sA.package_dir,'resources','mismatch_table_count.csv'),index_col=0)
table_mismatch_nucleotides = pd.read_csv(os.path.join(sA.package_dir,'resources','mismatch_table_mRNA_nucleotides.csv'),index_col=0).fillna('None')
table_mismatch_positions = pd.read_csv(os.path.join(sA.package_dir,'resources','mismatch_table_positions.csv'),index_col=0).fillna('None')
for i in table_mismatch_nucleotides.index:
    for c in table_mismatch_nucleotides.columns:
        table_mismatch_nucleotides.loc[i,c] = ast.literal_eval(table_mismatch_nucleotides.loc[i,c])
        table_mismatch_positions.loc[i,c] = ast.literal_eval(table_mismatch_positions.loc[i,c])
        
#Use lookup tables. Note that these tables assume L/I ambiguity, which are looked up under L. 
dfsubs['# BP Mismatches'] = dfsubs.apply(lambda row: table_mismatch_count.loc[row['Codon'],sA.dict_321ambiguous.get(row['Destination'])],axis=1)
dfsubs['BP Mismatch Nucleotides'] = dfsubs.apply(lambda row: table_mismatch_nucleotides.loc[row['Codon'],sA.dict_321ambiguous.get(row['Destination'])],axis=1)
dfsubs['# BP Mismatch Positions'] = dfsubs.apply(lambda row: table_mismatch_positions.loc[row['Codon'],sA.dict_321ambiguous.get(row['Destination'])],axis=1)

"""
Old way before the lookup table
progress_bar = tqdm(total=len(dfsubs.index))

dfsubs[['# BP Mismatches','BP Mismatch Positions',' BP Mismatch Nucleotides']] = dfsubs.apply(lambda row:
            lowest_bp_mismatch(row),axis=1,result_type='expand')
# Close the progress bar
progress_bar.close()
"""

#%%

sA.timestamp('Getting PSM rollup quantifications...')

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

protbySample['Sample'] = protbySample['File'].str.extract(r'(.*?)_')
protbySample['Replicate'] = protbySample['File'].str.extract(r'_(\d*)')

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
sampleQuant['Protein SD'] = bySample['Protein Intensity'].std()
sampleQuant['Protein CV'] =  sampleQuant['Protein SD']/sampleQuant['Protein Intensity']
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

