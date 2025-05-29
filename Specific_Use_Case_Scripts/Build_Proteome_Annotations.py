#A script to build lookup tables for protein values across the proteome. 
# Uses API for UniProt, AlphaFold to get positional annotations for each protein in a FASTA
# Uses DSSP for secondary structure and solvent accessibility
# Uses xxx for folding stability calculations [TBI]
# Uses yyy for evolutionary homology evaluation [Probably won't be implemented but would be fun]
# tjl 2024-10


#%% Variables to manually change
update = False  #Whether to lookup and overwrite values for proteins already in local database
    #I had to tell the code where to find the components.cif file
os.environ['LIBCIFPP_DATA_DIR'] = r'C:\Users\Champion Lab\source\repos\dssp\out\install\x64-Debug\share\libcifpp'
    #And where the mmcif dictionary was
mmcif_dict_path = r'C:\Users\Champion Lab\source\repos\dssp\out\install\x64-Debug\bin\mmcif_pdbx.dic'
#%% Package imports
import pandas as pd
import os
import re 
from tqdm import tqdm
from substitutionannotation.resources import assets as sA
from Bio import SeqIO,PDB
from Bio.PDB.DSSP import DSSP
from io import StringIO
import requests
import tempfile

#%% Existing lookup table, proteome import
path_structures = os.path.join(sA.package_dir,'resources','Protein Structure from Uniprot.csv')
path_solvent = os.path.join(sA.package_dir,'resources','Solvent Accessibility from DSSP.csv')
path_binding_site = os.path.join(sA.package_dir,'resources','Active or Binding Sites from Uniprot.csv')
path_signal = os.path.join(sA.package_dir,'resources','Signal Sequences from Uniprot.csv')
path_disulfide_bond = os.path.join(sA.package_dir,'resources','Disulfide bonds from Uniprot.csv')
path_af_structures = os.path.join(sA.package_dir,'resources','Structures from AlphaFold.csv')
try:
    structures = pd.read_csv(path_structures,low_memory=False)
    af_structures = pd.read_csv(path_af_structures,low_memory=False)
    solvent = pd.read_csv(path_solvent,low_memory=False)
    binding = pd.read_csv(path_binding_site,low_memory=False)
    signal = pd.read_csv(path_signal,low_memory=False)
    disulfide = pd.read_csv(path_disulfide_bond,low_memory=False)
except:
    raise Exception('Could not import local databases! Please check seed database with index of appropriate length.')

#Note I manually checked that the index length was longer than any protein sequence in my database.
#If you're adopting this code, you may need to account for that, as I manually assign attributes with the index
# representing the protein position.


#Get UniProt accessions from fasta
protein_index = SeqIO.index(sA.path_protein_fasta,format='fasta')
accessions = []
no_accessions = []      #In my fasta, there are 17 pseudogenes without a UniProt reference. 
                        #Most of these do have a Uniprot entry, but the sequences aren't 1:1
for key,rec in protein_index.items():
    accession = re.findall(r'UniProtKB/Swiss-Prot:(.+?)\]',rec.description)
    if accession:
        accession = accession[0]
        accessions.append(accession)
    else:
        no_accessions.append(key)


#Dictionary to convert one letter codes to verbose secondary structure
dict_dssp_readable = dict(H = 'α-helix',
B = 'residue in isolated β-bridge',
E = 'extended strand, participates in β ladder',
G = '310-helix',
I = 'π-helix',
P = 'κ-helix (poly-proline II helix)',
T = 'hydrogen-bonded turn',
S = 'bend')

#%% Loop through proteome, assign new values
progress_bar = tqdm(total=len(accessions))
for i,uniprot_accession in enumerate(accessions):
    progress_bar.update()
    if update == False:     #Check if we already have the UniProt ID in the local database
        if solvent.columns.str.contains(uniprot_accession).any():   
            continue
    
    #Pull Uniprot data
    response_uniprot = requests.get(f'https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_accession}')        
    if response_uniprot.status_code == 200:
        features = response_uniprot.json().get("features", [])
        for feature in features:
            try:
                cat = feature.get('category')
                begin = int(feature.get('begin'))
                end = int(feature.get('end'))
            except:         #Currently all the featuers I am interested in have these values
                pass            #I may have to adjust this if that changes
            if feature.get('type') == 'AlphaFoldDB':
                AF_ID = feature.get('id')
            if cat == 'STRUCTURAL':
                structures.loc[begin:end,uniprot_accession] = feature.get('type')
            if cat == 'DOMAINS_AND_SITES':
                if feature.get('type') in ['BINDING','ACT_SITE']:
                    binding.loc[begin:end,uniprot_accession] = True
                    pass
            if cat == 'PTM':
                if feature.get('type') == 'DISULFID':
                    disulfide.loc[begin,uniprot_accession] = True
                    disulfide.loc[end,uniprot_accession] = True
            if cat == 'MOLECULE_PROCESSING':
                if feature.get('type') == 'SIGNAL':
                    signal.loc[begin:end,uniprot_accession] = True

    else:
        raise Exception(f'Unable to get UniProt response for {uniprot_accession} \n Response code {response_uniprot.status_code}')
                        
    #There are two isoforms that may cause reference issues... I don't have protein evidence of them
    #so I'm opting to ignore the problem.
    if uniprot_accession == 'P63284-2':         #clpB has a truncated isoform, so we have to manually re-align positions from Uniprot
        pass
    if uniprot_accession == 'P02919-2':         #mrcB has a truncated isoform
        pass
    pass

    #AlphaFold structure or other PDB file can be used both to calculate secondary structure and other things
    response_alphafold = requests.get(f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}')
    if response_alphafold.status_code == 404:
        continue
    elif response_alphafold.status_code != 200:
        raise Exception(f'AlphaFold search failed! Status code {response_alphafold.status_code}')
    pdb_url = response_alphafold.json()[0].get('pdbUrl')
    response_pdb = requests.get(pdb_url)
    pdb_data = response_pdb.text
    parser = PDB.PDBParser()
    structure = parser.get_structure(uniprot_accession, StringIO(pdb_data))        # Extract the pLDDT scores (from the B-factor field)
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

        # Extract secondary structure,solvent accessibility information
            #See Bio.PDB.DSSP documentation for syntax explaination
        secondary_structure = [dssp[key][2] for key in dssp.keys()]
        relative_solvent = [dssp[key][3] for key in dssp.keys()]

    secondary_structure = pd.Series(secondary_structure)
    secondary_structure = secondary_structure.replace(dict_dssp_readable).replace('-','None')
    secondary_structure[mask_plddt] = 'NC'      #NC for No confidence, as opposed to confidence in non-structural segment
    relative_solvent = pd.Series(relative_solvent,dtype=object)
    relative_solvent[mask_plddt] = 'NC'
    os.remove(temp_pdb.name)
        #Assign values to proteome frame
    af_structures[uniprot_accession] = secondary_structure
    solvent[uniprot_accession] = relative_solvent
        
    #Assigning these values piecewises fragments the physical memory location of the frame, causing a performance warning
    #This resets the dataframe in memory space I guess?
    
    if i%100 == 0:
        structures = structures.copy()
        af_structures = af_structures.copy()
        solvent = solvent.copy()
        binding = binding.copy()
        disulfide = disulfide.copy()
        signal = signal.copy()
        domains = domains.copy()
progress_bar.close()
    
#%%  Write updated tables

structures.to_csv(path_structures)
af_structures.to_csv(path_af_structures)
solvent.to_csv(path_solvent)
binding.to_csv(path_binding_site)
signal.to_csv(path_signal)
disulfide.to_csv(path_disulfide_bond)
