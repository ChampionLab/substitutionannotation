#A script to quantify substitutions:
#A protein-centric quantification, summing all signals by different classifications for each aa position in a protein, for every protein
#Also reports some aggregate substitution metrics

#To be run AFTER filtering psms for known genomic mutations (such as Filter_Xac_genomic_mutations.py)
#Note: this uses regex that assumes no "_" in sample names. Otherwise, manual parsing of sample/replicate annotations is needed.
# Helpful protein ID for individual checks - lcl|U00096.3_prot_AAC76345.1_3255
#tjl 2024-06-14

#%% Package Imports
import pandas as pd
import numpy as np
import os
import re 
import datetime
from io import StringIO
from Bio import SeqIO,PDB
from Bio.PDB.DSSP import DSSP
from tqdm import tqdm
from substitutionannotation.resources import assets as sA
from Bio import SeqIO
import warnings
import math
import pickle
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from matplotlib.ticker import StrMethodFormatter
from fractions import Fraction
import datetime
import requests
import tempfile


#%%     Things to manually change
    #A regex tag to look for in a protein label, removing these proteins from the psm list
prot_exclude = r'_cont'

#%%
#Temporal uninteresting things
try:
    pd.set_option('future.no_silent_downcasting', True) #Silence a future warning from pandas v2.2.2
except:
    pass#Hurray you are using a more up to date version of pandas

#Dictionary to convert one letter codes to verbose secondary structure
dict_dssp_readable = dict(H = 'α-helix',
B = 'residue in isolated β-bridge',
E = 'extended strand, participates in β ladder',
G = '310-helix',
I = 'π-helix',
P = 'κ-helix (poly-proline II helix)',
T = 'hydrogen-bonded turn',
S = 'bend')

#%%     #Define a class to count things for each protein. My first usage of classes in code...

class CodonSpectralCounter:
    """A class for making/managing dfs that keep spectral counts by protein
    Initiate with an example column to set the level of multiindex (used to manage sample/replicate/aggregate
    distinctions)"""
    def __init__(self,columns):
        self.df = pd.DataFrame(columns = columns)                #DataFrame for storing stuff

    def set_dna_sequence(self,sequence):
        """Takes a dna sequence, splits it into codons, then assigns the column ['DNA Sequence'] with each row being one codon.
        Note: only compatible with 1-D set of columns!!"""
        codons = [sequence[i:i+3] for i in range(0,len(sequence),3)]
        self.df['DNA Sequence'] = codons

    def set_aa_sequence(self,sequence):
        """Takes an amino acid sequence, and assigns the column ['aa Sequence'] with one aa per row.
        Note: only compatible with 1-D set of columns!!"""
        aa = [aa for aa in sequence]
        self.df['aa Sequence'] = aa

    def add_count(self,counttype, countpos,n):
        """Each time this is called, adds n to the values of column counttype
        from row coutpos[0] to countpos[1] (inclusive).
        If counttype is new, creates the column, IF the new column matches the dimensions of the existing index/multiindex.
        
        """
        if not counttype in self.df.columns:        
            # Check if counttype dimensions match the multi-index columns of self.df
            expected_levels = self.df.columns.nlevels
            if type(counttype) == str:
                actual_levels = 1
            else:
                actual_levels = len(counttype)

            if expected_levels != actual_levels:
                raise ValueError(f"'{counttype}' dimensions do not match the column dimensions!")
        
            #Initialize new columns
            self.df[counttype] = float(0)      

        for bound in countpos:
            if not bound in self.df.index:               #Warn if adding counts out of sequence bounds
                warnings.warn('Tried to add count outside bounds of sequence!')

        self.df.loc[countpos[0]:countpos[1],counttype] = self.df.loc[countpos[0]:countpos[1],counttype] + n

    def get_count(self,counttype, position):
        return self.df.loc[position,counttype]

codons_by_aa_alphabetical = ['GCT', 'GCC', 'GCA', 'GCG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
                              'AAT', 'AAC', 'GAT', 'GAC', 'TGT', 'TGC', 'CAA', 'CAG', 'GAA', 'GAG',
                                'GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC', 'ATT', 'ATC', 'ATA', 'TTA',
                                  'TTG', 'CTT', 'CTC', 'CTA', 'CTG', 'AAA', 'AAG', 'ATG', 'TTT', 'TTC',
                                    'CCT', 'CCC', 'CCA', 'CCG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',
                                      'ACT', 'ACC', 'ACA', 'ACG', 'TGG', 'TAT', 'TAC', 'GTT', 'GTC', 'GTA',
                                        'GTG', 'TAA', 'TAG', 'TGA']
codons_by_aa_alphabetical_nostop = codons_by_aa_alphabetical[:-3]




codon_to_aa = {
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
    'AAT': 'Asn', 'AAC': 'Asn',
    'GAT': 'Asp', 'GAC': 'Asp',
    'TGT': 'Cys', 'TGC': 'Cys',
    'CAA': 'Gln', 'CAG': 'Gln',
    'GAA': 'Glu', 'GAG': 'Glu',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
    'CAT': 'His', 'CAC': 'His',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile',
    'TTA': 'Leu', 'TTG': 'Leu', 'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'AAA': 'Lys', 'AAG': 'Lys',
    'ATG': 'Met',
    'TTT': 'Phe', 'TTC': 'Phe',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser', 'AGT': 'Ser', 'AGC': 'Ser',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'TGG': 'Trp',
    'TAT': 'Tyr', 'TAC': 'Tyr',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'TAA': 'STOP', 'TAG': 'STOP', 'TGA': 'STOP'
}

# %% 

#Make output directory 
outputdir = os.path.join(sA.outputdir,'Substitution Quant Output')
if not os.path.exists(outputdir):
    os.mkdir(outputdir)

#Import data
print('Importing data')
print(datetime.datetime.now())
progress_bar = tqdm(desc='Fake loading bar!')
progress_bar.close()
allpsms = pd.read_csv(os.path.join(sA.outputdir,'AllPSMsAndFilters.csv')).infer_objects()
#Reassign solely numerical Sample IDs 
allpsms['Sample'] = allpsms['Sample'].astype(str)
samples = allpsms['Sample'].unique()
#Check if there were replicates
if allpsms['Replicate'].isna().all():
    allpsms['Replicate'] = allpsms['Replicate'].fillna(1).astype(int)
    
sample_reps = (allpsms['Sample'] + '_' + allpsms['Replicate'].astype(str)).unique()
fragger_protein = pd.read_csv(os.path.join(sA.inputdir,"combined_protein.tsv"),sep='\t').infer_objects()

#Check if known genomic mutations were filtered upstream
try:
    ionquant = pd.read_csv(os.path.join(sA.outputdir,'Ionquant Peptide.csv'))
except:
    ionquant = pd.read_csv(os.path.join(sA.inputdir,'combined_modified_peptide.tsv'),sep='\t',low_memory=False).infer_objects()

#Fill nan intensities with 0's
allpsms['Intensity'] = allpsms['Intensity'].fillna(0)
allpsms = allpsms.drop_duplicates(subset='Spectrum')    #allpsms has duplicated spectra for Leu/Ile destination substitutions, we don't want to double count these
                                                        #Note this means you should not be using 'Substituted Sequence' for alignments, as you lose one of the possibilities
#Filter out proteins
allpsms = allpsms[~allpsms['Protein'].str.contains(prot_exclude)]


#Parse substitution PSMs
subpsms = allpsms[allpsms['Is Sub']]
subpsms = subpsms[~subpsms['Is Danger']]

#Get substitution position
def get_substitution_position(row):
    if row['Is Sub']:
        pos = re.search(r'.\[',row['Modified Peptide']).start() + row['Protein Start']
        return pos
    else:
        return None

allpsms['Substitution Position'] = allpsms.apply(get_substitution_position,axis=1)

#%%
print('Getting uniprot IDs')
#Get protein:uniprot dictionary from search
dict_uniprot = {}
unique_prot = fragger_protein.drop_duplicates(subset='Protein')
unique_prot['Uniprot ID'] = unique_prot['Description'].str.extract(r'\[db_xref=UniProtKB/Swiss-Prot:(\w+)\]')
#Look up protein:uniprot IDs from manual table
df_uniprot_map = pd.read_csv(os.path.join(sA.package_dir,r'resources\ManualUniProt ID map.csv'))
df_uniprot_map = df_uniprot_map.set_index('GenBank ID')
for id in unique_prot.loc[unique_prot['Uniprot ID'].isna(),'Protein']:
    if id in df_uniprot_map.index:
        unique_prot.loc[unique_prot['Protein'] == id, 'Uniprot ID'] = df_uniprot_map.loc[id,'UniProt ID']

n_missing = unique_prot['Uniprot ID'].isna().sum()
print(f'{n_missing} proteins did not have a mapped Uniprot ID!')
print('These IDs are stored in the no_uniprot variable, but these proteins will not be quantified...')
no_uniprot = unique_prot.loc[unique_prot['Uniprot ID'].isna(),'Protein']
unique_prot = unique_prot.dropna(subset='Uniprot ID')

#Make a dictionary to look up uniprot Ids from protein IDs
for _,row in unique_prot.iterrows():
    dict_uniprot[row['Protein']] = row['Uniprot ID']

#Parse peptide and substitution LFQ
    #A couple of aligning names/annotations
        #And replacement for files with no replicates
dict_noreplicate = {}
for sample in samples:
    dict_noreplicate[sample] = sample+'_1'
    
    #Strategy is to extract columns with exact Sample_Replicate names - strip the extra annotations from the target columns
ionquant.columns = ionquant.columns.str.replace(' MaxLFQ Intensity','')     #Strip MaxLFQ to extract
ionquant = ionquant.rename(columns=dict_noreplicate)
ionquant_type = ionquant.copy()
ionquant_type.columns = ionquant_type.columns.str.replace(' Match Type','')     #Strip Match Type to extract
ionquant = ionquant.rename(columns={'Modified Sequence':'Modified Peptide'})
ionquant_type = ionquant_type.rename(columns={'Modified Sequence':'Modified Peptide'})
ionquant.loc[:,'Modified Peptide'] = ionquant['Modified Peptide'].str.replace(r'C\[57.0215\]','',regex=True) #Remove fixed carbamidomethylation annotations
ionquant_type.loc[:,'Modified Peptide'] = ionquant_type['Modified Peptide'].str.replace(r'C\[57.0215\]','',regex=True) #Remove fixed carbamidomethylation annotations

    #Reshape data to long format
        #Direct Quant/ MBR long format
ionquant_type = ionquant_type.melt(id_vars=['Modified Peptide','Protein','Start','End','Protein ID'],value_vars= sample_reps,value_name='Match Type',var_name='Sample_Rep')
ionquant_type = ionquant_type[~(ionquant_type['Match Type'] == 'unmatched')]  #Drop 0 value assignments
ionquant_type.loc[:,'Replicate'] = ionquant_type['Sample_Rep'].str.extract(r'_(.*)').astype(int)
ionquant_type.loc[:,'Sample'] = ionquant_type['Sample_Rep'].str.extract(r'(.*)_')
ionquant_type['Is MBR'] = ionquant_type['Match Type'] == 'MBR'
ionquant_type = ionquant_type[['Sample','Replicate','Modified Peptide','Is MBR']]
        #All quant data, long format, and renaming columns
allLFQ = ionquant.melt(id_vars=['Modified Peptide','Protein','Start','End','Protein ID'],value_vars= sample_reps,value_name='LFQ',var_name='Sample_Rep')
allLFQ['Gene'] = allLFQ['Protein ID'].str.extract(r'gene=([^]]+)')                              #Get gene name
allLFQ = allLFQ.loc[allLFQ['LFQ'] > 0,:]                                                        #Drop 0 signal peptide/sample/rep combos
allLFQ.loc[:,'Replicate'] = allLFQ['Sample_Rep'].str.extract(r'_(.*)').astype(int)
allLFQ.loc[:,'Sample'] = allLFQ['Sample_Rep'].str.extract(r'(.*)_')
allLFQ = allLFQ.drop(columns = ['Sample_Rep'])
allLFQ = allLFQ.rename(columns={'Start':'Protein Start','End':'Protein End'})
allLFQ = allLFQ.merge(ionquant_type,how='left',on=['Sample','Replicate','Modified Peptide'])    #Merge in MBR info
allLFQ['Is Sub'] = allLFQ['Modified Peptide'].isin(subpsms['Modified Peptide'])                 #Align which quantified sequences represent substitutions
   
df_ann_rep = allpsms.sort_values(by='Substitution Position').drop_duplicates(subset=['Modified Peptide'])
allLFQ = allLFQ.merge(df_ann_rep[['Modified Peptide','Destination','Substitution Position']],how='left',on=['Modified Peptide'])


#Get Sampe/Replicate/Quant index combos for initializing dfs
replicates = allpsms['Replicate'].unique()
sample_rep_sets = allpsms.groupby(by=['Sample','Replicate'])['Is Sub'].count().index
sample_rep_quant_idx_values = []
# Make 3rd tier to avoid non-existant Sample-Rep pairs
for combo in sample_rep_sets:
    sample_rep_quant_idx_values.extend([(combo[0], combo[1],'Total')])
    sample_rep_quant_idx_values.extend([(combo[0], combo[1],'Any Substitution')])

# Create a new multi-index with the new level
sample_rep_quant_cols = pd.MultiIndex.from_tuples(sample_rep_quant_idx_values, names=sample_rep_sets.names + ['Quant'])

#Filter MBR quantifications     --I hope this is temporary, until I can get reasonable results from MBR--
mbrLFQ = allLFQ[allLFQ['Is MBR']]
allLFQ = allLFQ[~allLFQ['Is MBR']]                                      #Remove ALL MBR quantifications

#Write tables
print('Writing peptide-centric quantification tables')
print(datetime.datetime.now())

allLFQ.to_csv(os.path.join(sA.outputdir,'Allpep IonQuant.csv'))
mbrLFQ.to_csv(os.path.join(sA.outputdir,'Removed MBR quantifications by IonQuant.csv'))


#%%  
#Calculate the protein-centric quantifications
print('Getting protein-centric quantifications')
print(datetime.datetime.now())

#Import genome, set up variables to store values
print('Initializing protein dataframes will be slower if script needs to find secondary strucutre information from UniProt')
print('This will be saved to a file to speed up future runs.')
genome_index = SeqIO.index(sA.path_dna_fasta,format='fasta')
protein_index = SeqIO.index(sA.path_protein_fasta,format='fasta')
structures = pd.read_csv(os.path.join(sA.package_dir,'resources','Structures from AlphaFold.csv'),low_memory=False)

#Initialize a dictionary with the accession ID as the key and a CodonSpectralCounter item as the value
dict_area_sample_rep_proteins = {}   #No aggregation by sample/replicate
#Variables to store annotation failures
no_uniprot = []
no_structure = []
#Get protein sequences in aggregate dfs
progress_bar.close()
progress_bar = tqdm(total=len(unique_prot['Protein']),desc='Initializing protein dataframes')
for prot in unique_prot['Protein']:
    progress_bar.update()
    has_uniprot = True
    if dict_area_sample_rep_proteins.get(prot):         #Test for duplicate entries
                    warnings.warn(f'Duplicate {prot} entries. Previous entry will be overwritten')
    if not prot in genome_index.keys():
        warnings.warn(f'Identified protein {prot} was not in the reference FASTA! Is this a contaminant?')
        continue
    
        #Initialize Total columns
    counter_area = CodonSpectralCounter(sample_rep_quant_cols)
        #Get index information
    gene = prot.replace('_prot_','_cds_')
    seq = genome_index.get(gene).seq
    aa = protein_index.get(prot).seq
    dnaseq = str(seq)              
            #Check that DNA/protein lengths align
    if len(aa)* 3 != len(dnaseq):
        if dnaseq[-3:] in ['TGA','TAG','TAA']:
            dnaseq = dnaseq[:-3]              #Trim off stop codon when necessary
            if len(aa)* 3 != len(dnaseq):
                warnings.warn(f'Mismatch between DNA and protein sequence length for {gene}')
    dnalen = len(dnaseq)
    if dnalen % 3 != 0:                       #Test if gene is proper length for translation
        warnings.warn(f'{prot} had a partial codon! This will be trimmed') 
        end = math.floor(dnalen / 3) * 3
        dnaseq = str(seq)[0:end]

    position = range(0,len(aa))
    aa = [letter for letter in aa]
    codon = [dnaseq[i:i+3] for i in range(0,len(dnaseq),3)]
    relpos = pd.cut(position, bins=100, include_lowest=True, labels=np.linspace(0, 1, 100))
        #Get secondary structure information
    try:        #Get uniprot from dictionary
        uniprot = dict_uniprot.get(prot)
    except KeyError:            #Get uniprot from FASTA
        uniprot = re.findall(r'UniProtKB/Swiss-Prot:(.+?)\]',protein_index.get(prot).description)[0]
    if not uniprot:
            #Get uniprot from FASTA, not sure why this sometimes returns none instead of key error
        try:
            uniprot = re.findall(r'UniProtKB/Swiss-Prot:(.+?)\]',protein_index.get(prot).description)[0]
        except IndexError:
            no_uniprot.append(prot)
            has_uniprot = False

    if has_uniprot:
        try:        #Get structure from pre-compiled dataframe of structures
            secondary_structure = structures[uniprot]
        except KeyError:
            #AlphaFold structure or other PDB file can be used both to calculate secondary structure and other things
            response_alphafold = requests.get(f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot}')
            if response_alphafold.status_code == 404:       #Ignore proteins with no prediction by alphafold
                
                secondary_structure = pd.Series(index=range(len(aa)),data = ['NA']*len(aa))        #Todo: Run local alphafold on these proteins? They are only ones that are very short
            elif response_alphafold.status_code != 200:     
                no_structure.append(uniprot)
                secondary_structure = pd.Series(index=range(len(aa)),data = ['NA']*len(aa))
                
            else:           #Use DSSP to get the secondary structure from the AlphaFold PDB file
                pdb_url = response_alphafold.json()[0].get('pdbUrl')
                response_pdb = requests.get(pdb_url)
                pdb_data = response_pdb.text
                parser = PDB.PDBParser()
                structure = parser.get_structure(uniprot, StringIO(pdb_data))        # Extract the pLDDT scores (from the B-factor field)
                plddt_scores = []
                for model in structure:  # Iterate through the model(s)
                        for chain in model:  # Iterate through chains
                            for residue in chain:  # Iterate through residues (amino acids)
                                # Check if the residue has a Cα (CA) atom
                                if 'CA' in residue:
                                    # Append the B-factor (pLDDT score) of the Cα atom
                                    plddt_scores.append(residue['CA'].bfactor)
                plddt_scores = pd.Series(plddt_scores)
                mask_plddt = plddt_scores < 70  #A mask for 'confident' AF structures. Can adjust this cutoff as desired.
                    #Use a temporary file to store the pdb and calculate structure stuff
                with tempfile.NamedTemporaryFile(suffix=".pdb",delete=False) as temp_pdb:
                    temp_pdb.write(pdb_data.encode())  # Write PDB data as bytes
                    temp_pdb.flush()  # Ensure the data is written to the file

                    # Use DSSP on the temporary file
                    model = structure[0]  # First model in the structure
                    dssp = DSSP(model, temp_pdb.name,'mkdssp')  # Use the file path

                    # Extract secondary structure
                        #See Bio.PDB.DSSP documentation for syntax explaination
                    secondary_structure = [dssp[key][2] for key in dssp.keys()]

                secondary_structure = pd.Series(secondary_structure)
                secondary_structure = secondary_structure.replace(dict_dssp_readable).replace('-','None')
                secondary_structure[mask_plddt] = 'NC'      #NC for No confidence, as opposed to confidence in non-structural segment
                os.remove(temp_pdb.name)
                    #Assign values to proteome frame
                structures[uniprot] = secondary_structure
             

        #stored structural data has a vector longer than all proteins. Trim excess
    structure = secondary_structure[:len(aa)]
        #Newly pulled structure information does not always reach the C-terminus. Append full length annotation
    tail = pd.Series(index=range(len(structure),len(aa)))
    structure = pd.concat([structure,tail]).fillna('NA')

    df_index = pd.MultiIndex.from_arrays([position,codon,aa,relpos,structure],names=['Position','Codon','AA','Relative Position','Structure'])
    df_init = pd.DataFrame(index=df_index,columns=sample_rep_quant_cols).fillna(0.0)
    counter_area.df = df_init.copy()
    dict_area_sample_rep_proteins[prot] = counter_area

progress_bar.close()
if no_uniprot:
    print(f'Error in getting uniprot ID for {len(no_uniprot)} proteins')
if no_structure:
    print(f'Error in getting structure for {len(no_structure)} proteins')

#Update precompiled structures table
structures.to_csv(os.path.join(sA.package_dir,'resources','Structures from AlphaFold.csv'),index=False)

#%%     Sum LFQ values by protein
lfq_by_prot = allLFQ.groupby(by='Protein')

def sum_lfq(df_recipient, row):
    df_recipient.loc[row['Protein Start']-1:row['Protein End']-1,(row['Sample'],row['Replicate'],'Total')] += row['LFQ']
    if row['Is Sub']:
        df_recipient.loc[row['Substitution Position']-1,(row['Sample'],row['Replicate'],'Any Substitution')] +=row['LFQ']
        df_recipient.loc[row['Substitution Position']-1,(row['Sample'],row['Replicate'],row['Destination'])] += row['LFQ']

progress_bar = tqdm(total=len(lfq_by_prot.groups.keys()))
progress_bar.set_description(f'Summing LFQ values by protein...')
for prot in lfq_by_prot.groups.keys():
    progress_bar.update()
    if dict_area_sample_rep_proteins.get(prot) == None:
        warnings.warn(f'No protein {prot} initialized!')
        continue
        
    dff = lfq_by_prot.get_group(prot)
    len_prot = len(protein_index.get(prot).seq)
    unique_combinations = dff[['Sample', 'Replicate', 'Destination']].dropna(subset='Destination').drop_duplicates()
    destinations_columns = list(unique_combinations.itertuples(index=False, name=None))
    merge_cols = pd.MultiIndex.from_tuples(destinations_columns + sample_rep_quant_idx_values,names = ('Sample','Replicate','Quant'))  
    df_merge = pd.DataFrame(index=range(0,len_prot),columns = merge_cols).fillna(0)
    dff.apply(lambda x: sum_lfq(df_merge,x), axis =1)
    df_merge.index = dict_area_sample_rep_proteins.get(prot).df.index
    dict_area_sample_rep_proteins.get(prot).df = df_merge
progress_bar.close()

#%%         Export data
print('Exporting individual protein quant info')
print(datetime.datetime.now())
outputdir = os.path.join(sA.outputdir,'Substitution Quant Output')
if not os.path.exists(outputdir):
    os.mkdir(outputdir)

with open(os.path.join(outputdir,'dict_protein_sample_rep_counters_area.pkl'), 'wb') as file:
    pickle.dump(dict_area_sample_rep_proteins, file)

# %%   
# Get whole proteome table
print('Protein-oriented quant complete! Combining protein tables...')
print(datetime.datetime.now())


    #Concatenate all protein quant
point_rep_area = pd.concat([entry.df for entry in dict_area_sample_rep_proteins.values()],keys = dict_area_sample_rep_proteins.keys(),names =['Protein'])
point_rep_area = point_rep_area.sort_index().sort_index(axis=1).replace(0,np.nan)

    #Get normalized values
print('Converting to long format...')
print(datetime.datetime.now())
df_rep_long = point_rep_area.stack(level=['Sample','Replicate','Quant']).rename('Signal').reset_index()
df_long_totals = df_rep_long.loc[df_rep_long['Quant'] == 'Total'].drop(columns='Quant')
df_rep_long = df_rep_long.loc[~(df_rep_long['Quant'] == 'Total'),['Protein','Position','Sample','Replicate','Signal','Quant']]
df_long_totals.rename(columns={'Signal':'Total'},inplace=True)
df_long = df_rep_long.merge(df_long_totals,how='outer',on=['Protein','Position','Sample','Replicate'])
df_long.rename(columns={'Signal':'Substitution Signal'},inplace=True)




df_long['Total'] = np.log10(df_long['Total'].astype(float))
df_long['Substitution Signal'] = np.log10(df_long['Substitution Signal'].astype(float))
df_long['Substitution Proportion'] = df_long['Substitution Signal'] - df_long['Total']
df_long['Substitution %'] = 10 ** df_long['Substitution Proportion'] * 100

point_rep_area_norm = df_long.pivot_table(index=['Protein','Position'],columns=['Sample','Replicate','Quant'],values='Substitution Proportion')

#Aggregate to location and incorporated aa
n_replicates = len(df_long['Sample'].unique()) * len(df_long['Replicate'].unique())
df_sum = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','AA'])[['Substitution Signal','Total','Substitution Proportion']].sum()
df_count = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','AA'])[['Substitution Signal','Total','Substitution Proportion']].count()
df_agg = df_sum/df_count
df_cvs = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','AA'])[['Substitution Signal','Total','Substitution Proportion']].std()/df_agg
df_cvs['Substitution Proportion'] = np.abs(df_cvs['Substitution Proportion'])
df_cvs.columns = df_cvs.columns + ' CV'
df_agg = df_agg.join(df_cvs)
df_agg['# Observations'] = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','AA'])['Substitution Signal'].count()
df_agg['Substitution %'] = 10 ** df_agg['Substitution Proportion'] * 100
df_agg = df_agg.reset_index()
df_agg.rename(columns={'AA':'Encoded aa','Quant':'Incorporated aa'},inplace=True)

#Aggregate to location, sample, and incorporated aa
n_samples = len(df_long['Sample'].unique())
df_sum = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','Sample','AA'])[['Substitution Signal','Total','Substitution Proportion']].sum()
df_count = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','Sample','AA'])[['Substitution Signal','Total','Substitution Proportion']].count()
df_sample = df_sum/df_count
df_cv_sample = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','Sample','AA'])[['Substitution Signal','Total','Substitution Proportion']].std()/df_sample
df_cv_sample['Substitution Proportion'] = np.abs(df_cv_sample['Substitution Proportion'])
df_cv_sample.columns = df_cv_sample.columns + ' CV'
df_sample = df_sample.join(df_cv_sample)
df_sample['# Observations'] = df_long[~(df_long['Quant']=='Any Substitution')].groupby(by=['Protein','Position','Quant','Sample','AA'])['Substitution Signal'].count()
df_sample['Substitution %'] = 10 ** df_sample['Substitution Proportion'] * 100
df_sample = df_sample.reset_index()
df_sample.rename(columns={'AA':'Encoded aa','Quant':'Incorporated aa'},inplace=True)


#%%
print('Saving normalized and aggregated tables...')
if not os.path.isdir(os.path.join(outputdir,'Aggregated Tables')):
    os.mkdir(os.path.join(outputdir,'Aggregated Tables'))
   
point_rep_area.to_csv(os.path.join(outputdir,'Aggregated Tables','point_rep_area.csv'))
point_rep_area_norm.to_csv(os.path.join(outputdir,'Aggregated Tables','point_rep_area_norm.csv')) 
df_long.to_csv(os.path.join(outputdir,'Aggregated Tables','df_rep_long.csv'))
df_agg.to_csv(os.path.join(outputdir,'Aggregated Tables','df_agg.csv'))
df_sample.to_csv(os.path.join(outputdir,'Aggregated Tables','df_sample.csv'))



#%%
#Plot an individual protein
#EfTu
prot = 'lcl|U00096.3_prot_AAC76364.1_3274'
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
#prot = 'sp|WOO|Dual_Luciferase'

indv_prot = df_rep_long[df_rep_long['Protein'] == prot]
indv_prot = indv_prot.reset_index()
indv_prot.to_csv(os.path.join(outputdir,prot.replace('|',' ') + 'long substitutions.csv'))

indv_prot_bysub = indv_prot[~indv_prot['Quant'].str.contains('Any')]
#indv_prot_bysub['Substitution Frequency'] = 10**indv_prot_bysub['Substitution Frequency'] 
indv_prot_anysub = indv_prot[indv_prot['Quant'].str.contains('Any')]
#indv_prot_anysub['Substitution Frequency'] = 10**indv_prot_anysub['Substitution Frequency'] 

indv_prot_total = point_sample_area.xs(prot,axis=0,level='Protein').xs('Total',axis=1,level='Quant')
indv_prot_total = indv_prot_total.stack()
indv_prot_total.name = 'Total Signal'
indv_prot_total = indv_prot_total.reset_index()

protlen = len(protein_index.get(prot).seq)
protseq = protein_index.get(prot).seq

dict_sample_colors = {'A':'red','X':'grey','E':'blue'}

fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.1)
#Substitution type subplot
for trace in px.scatter(
        indv_prot_bysub,
        x='Position',
        y='Substitution Frequency',
        color='Quant',
        labels={'Position': 'Position', 'Substitution Frequency': 'Substitution Frequency'}
    ).data:
    fig.add_trace(trace, row=1, col=1)
    
#Any substitution by sample subplot
for trace in px.scatter(
        indv_prot_anysub,
        x='Position',
        y='Substitution Frequency',
        color='Sample',
        color_discrete_sequence = ['red','blue','grey'],
        labels={'Position': 'Position', 'Substitution Frequency': 'Substitution Frequency'}
    ).data:
    fig.add_trace(trace, row=2, col=1)

#Total signal plot
for trace in px.scatter(
        indv_prot_total,
        x='Position',
        y='Total Signal',
        color='Sample',
        color_discrete_sequence = ['red','blue','grey'],
        labels={'Position': 'Position', 'Substitution Frequency': 'Substitution Frequency'}
    ).data:
    fig.add_trace(trace, row=3, col=1)


#log y axis
fig.update_yaxes(type="log", row=1, col=1, title_text='Indv Sub')
fig.update_yaxes(type="log", row=2, col=1, title_text='Any Sub')
fig.update_yaxes(type="log", row=3, col=1, title_text='Total Signal')
# Update x-axis to include custom ticks and labels
fig.update_xaxes(
    tickmode='array',
    tickvals=list(range(0,protlen)),
    ticktext=[a for a in protseq],
    range=[0, protlen],
    tickangle=0
)

for i in range(0, protlen, 25):
    fig.add_annotation(
        x=i,
        y=-0.1,  # Position well below the plot
        text=str(i),  # The number to display
        showarrow=False,
        yshift=0,  # No need for additional shifting
        xanchor='center',
        xref="x",
        yref="paper",  # Reference the paper coordinates for y (0 to 1 scale)
        font=dict(size=10)
    )
fig.update_layout(
    legend=dict(
        orientation="h",  # Horizontal orientation
        yanchor="bottom",  # Vertical anchor (bottom of the legend)
        y=1.02,  # Position above the plot area
        xanchor="center",  # Horizontal anchor (centered)
        x=0.5,  # Center the legend horizontally
        title='',  # Title for the legend (optional)
        traceorder="normal"  # Order of legend items
    ),
    legend_title_font=dict(size=12),  # Font size for legend title
    legend_font=dict(size=10)  # Font size for legend items
)
highlight = False
if highlight:
    #Specific box to highlight linker sequence
    fig.add_shape(
        go.layout.Shape(
            type="rect",
            x0=312, x1=329,    # Set x-range for the shaded region
            y0=10000, y1=10**8,        # Span the entire y-axis (relative range)
            xref="x3", yref="y3",   # xref = x-axis, yref = relative to paper
            fillcolor="LightGrey",  # Set fill color
            opacity=0.3,             # Set opacity of the fill
            line_width=0),              # No border line
        row=3, col =1
    )

fig.show(renderer='browser')
pio.write_html(fig,file=os.path.join(outputdir,prot.replace('|',' ') + 'long substitutions.html'))

