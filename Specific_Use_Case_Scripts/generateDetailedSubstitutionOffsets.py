#A script to write the detailed delta-mass for labile search in MSFragger
#2024-02-06 tjl

import pandas as pd
from substitutionannotation.resources import assets as sA
import os

detailed_substitution_offsets = pd.DataFrame(columns = ['# mass','Allowed sites','Diagnostic ions','Peptide remainder ions','Fragment remainder ions'],data=[[0,'','','','']])

#Melt delta mass matrix into long format
subs = sA.dfdm.reset_index().melt(id_vars='index',var_name='Origin',value_name='dm')
#Change 3 letters to 1 letter, align column names
subs['Allowed sites'] = subs['Origin'].replace(sA.dict_321)
#Remove identities from matrix
subs = subs[~(subs['dm'] == 0)]
#Expand I/L ambiguity
il = subs[subs['Origin'].str.contains('Ile')].reset_index()
subs = subs[~subs['Origin'].str.contains('Ile')]
il['Allowed sites'] = 'L'
subs = pd.concat([subs,il])
il['Allowed sites'] = 'I'
subs = pd.concat([subs,il])
#Assign dm to # mass [for precursor matching of mod in search]
subs['# mass'] = subs['dm']
#Assign dm to Fragment remainder ions [for MS2 matching and localization of mod in search]
subs['Fragment remainder ions'] = subs['dm']
#Empty columns
subs[['Diagnostic ions','Peptide remainder ions']] = ''


#Add in dangermods to list
danger = sA.dfdangermods.copy()
#Assign dm to # mass [for precursor matching of mod in search]
danger['# mass'] = danger['Mass Shift']
#Assign allowed sites, leaving blank any that are ambiguous on any residue
    #Note X denotes any reside in the danger reference, and has been generously applied to some PTMs that could be restricted
danger['Allowed sites'] = danger['Modified Residues'].str.replace(',','')
danger.loc[danger['Allowed sites'].str.contains('X'),'Allowed sites'] = ''
danger[['Diagnostic ions','Peptide remainder ions','Fragment remainder ions']] = '','',''

#Manually add in some fragment info for phospho modifications. Manual addition of other modification fragment details would also likely be beneficial
danger.loc[27,'Fragment remainder ions'] = -18.01056
danger.loc[27,'Allowed sites'] = 'ST'
danger.loc[30,:] =['Phosphorylation','Phosphorylation',79.966331,'STY',None,79.966331,'Y',216.043,'','']

detailed_substitution_offsets = pd.concat([detailed_substitution_offsets,subs,danger],join='inner')

detailed_substitution_offsets.to_csv(os.path.join(sA.package_dir,'resources','DetailedSubstitutionOffsets.tsv'),sep='\t',index=False)