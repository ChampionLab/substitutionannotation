# -*- coding: utf-8 -*-
"""
Created on Oct 23 2023
A script to create many substitutions distribution plots
Must follow FindSubs_xxx.py
@author: taylo
"""

#%%
#Import packages
import pandas as pd
import numpy as np
import seaborn as sns
import os
from matplotlib import pyplot as plt
from substitutionannotation.resources import assets as sA
import upsetplot
import re
import scipy.stats as st
import statsmodels.stats.multicomp as mc

#%%
#Import filtered data, add annotation data to quant data
dffiltered = pd.read_csv(os.path.join(sA.outputdir,'Filtered SSP Quant.csv'))
dfann = pd.read_csv(os.path.join(sA.outputdir,'SSP PSM.csv'))
anncols = ['# BP Mismatches','BP Mismatch Positions','Origin','Destination','PTMs','Codon']
ann = dfann.groupby(by=['Sample','Modified Peptide'])[anncols].min()
df = dffiltered.merge(ann.reset_index(),on=['Sample','Modified Peptide'],how='left')
df['logIntensity Base'] = np.log10(df['Base Intensity'].dropna())

dist_path = os.path.join(sA.outputdir,'Distribution Figures')
if not os.path.exists(dist_path):
    os.mkdir(dist_path)

# %%
#All sample intensity heatmap

hexplot = sns.jointplot(data=df,x='logIntensity Base',y='logIntensity',
                        kind='hex',gridsize=25,marginal_kws=dict(bins=25))
plt.subplots_adjust(left=0.2, right=0.8, top=0.8, bottom=0.2)  # shrink fig so cbar is visible
plt.ylabel('Substitution Intensity (log10)')
plt.xlabel('Genomic Cognate Intensity (log10)')
# make new ax object for the cbar
cbar_ax = hexplot.fig.add_axes([.85, .25, .05, .4])  # x, y, width, height
plt.colorbar(cax=cbar_ax,label='# of Substitutions')
plt.savefig(os.path.join(dist_path,'Filtered Intensity distributions'))
#%%
#Substitution Frequency by sample

#Find out if we need significance bars
data = df.pivot_table(columns='Sample',index='Modified Peptide',values='logRatio')
anova = [data[sample].dropna().to_list() for sample in data.columns]

anovaP = st.f_oneway(*anova)[1]

sns.violinplot(data = df, x = 'Sample', y = 'logRatio', cut = 0, density_norm='count')
plt.ylabel('Substitution Frequency')

if anovaP >0.05:
    plt.annotate('ANOVA P = {0:.3f}'.format(anovaP),(2.4,1.1))

plt.savefig(os.path.join(dist_path,'Frequency by Sample'))


# %%
#Substitution Frequency by codon


fig,ax = plt.subplots()
sns.violinplot(data = df, x = 'Codon', y = 'logRatio', cut = 0, density_norm='count')
plt.ylabel('Substitution Frequency')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

plt.savefig(os.path.join(dist_path,'Frequency by codon'))


# %%
#Substitution Frequency by destination


fig,ax = plt.subplots()
sns.violinplot(data = df, x = 'Destination', y = 'logRatio', cut = 0, density_norm='count')
plt.ylabel('Substitution Frequency')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.savefig(os.path.join(dist_path,'Frequency by destination'))


# %%
#Substitution frequency by origin
fig,ax = plt.subplots()
sns.violinplot(data = df, x = 'Origin', y = 'logRatio', cut = 0, density_norm='count')
plt.ylabel('Substitution Frequency')
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
plt.savefig(os.path.join(dist_path,'Frequency by origin'))


# %%
#Substitution mean Frequency by origin, destination

by = df.groupby(by=['Origin','Destination'])['logRatio'].mean().reset_index().pivot_table(columns='Destination',index='Origin',values='logRatio')
byann = df.groupby(by=['Origin','Destination'])['logRatio'].count().reset_index().pivot_table(columns='Destination',index='Origin',values='logRatio')

ax = sns.heatmap(data = by,cmap='plasma_r',cbar_kws={'label': 'log10 Substitution Frequency'}, annot=byann, square=True, fmt='.0f')
ax.set_facecolor((.2,.2,.2,.15))
plt.savefig(os.path.join(dist_path,'Frequency mean by Origin and Destination'))


# %%
#Substitution Frequency by # of bp mismatches

#Find out if we need significance bars
data = df.pivot_table(columns='# BP Mismatches',index='Modified Peptide',values='logRatio')
anova = [data[sample].dropna().to_list() for sample in data.columns]

anovaP = st.f_oneway(*anova)[1]

sns.violinplot(data = df, x = '# BP Mismatches', y = 'logRatio', cut = 0, density_norm='count')
plt.ylabel('Substitution Frequency')

if anovaP >0.05:
    plt.text(2.4,1.1,'ANOVA P = {0:.3f}'.format(anovaP), ha='right')

plt.savefig(os.path.join(dist_path,'Frequency by # BP Mismatches'))


# %%
#Substitution Frequency by bp mismatch position
#This is a hack at collapsing mismatch positions
df['Mismatch Positions'] = df['BP Mismatch Positions'].str.extract(r'([\d]{1,3})')
#Find out if we need significance bars
data = df.pivot_table(columns='Mismatch Positions',index='Modified Peptide',values='logRatio')
anova = [data[sample].dropna().to_list() for sample in data.columns]

anovaP = st.f_oneway(*anova)[1]

if anovaP < 0.05:
    tukey = data.reset_index().melt(id_vars='Modified Peptide',var_name='Mismatch Positions',value_name='logRatio').dropna(subset='logRatio')
    comp = mc.MultiComparison(data=tukey['logRatio'],groups = tukey['Mismatch Positions'])
    result = comp.tukeyhsd()
    dfresult = pd.DataFrame(data=result._results_table.data, columns=result._results_table.data[0]) #Store tukey results in DF because I can parse those
    dfresult = dfresult[dfresult['reject']==True]  #Pull out those with p < 0.05

fig,ax = plt.subplots()
sns.violinplot(data = df, x = 'Mismatch Positions', y = 'logRatio', cut = 0, density_norm='count', ax=ax)
plt.ylabel('Substitution Frequency')

if anovaP >0.05:
   plt.text(2.4,1.1,'ANOVA P = {0:.3f}'.format(anovaP), ha='right')

if anovaP < 0.05:
    if len(dfresult) > 0:
        #Get position of tick labels to line up significance indicators
        dict_ticks = {}
        for tick in ax.get_xticklabels():
            dict_ticks[tick.get_text()] = tick.get_position()[0]
        dfresult['group1'].replace(dict_ticks,inplace=True)
        dfresult['group2'].replace(dict_ticks,inplace=True)
        
        offset = 1.05 #Offset for significance lines
        lineh = df['logRatio'].max()
        for _,row in dfresult.iterrows():
            plt.hlines(lineh*offset,int(row['group1']),int(row['group2']), color='k')
            offset += 0.05

plt.savefig(os.path.join(dist_path,'Frequency by Mismatch Position'))



# %%
#Upset plot by sample


subsbysample = df.groupby(by='Sample')['Modified Peptide'].apply(set)

subsUp = upsetplot.from_contents({sample:subsbysample[sample] for sample in dfann['Sample']})

ax_dict = upsetplot.UpSet(subsUp, subset_size='count', min_subset_size=10,sort_by='cardinality').plot()
ax_dict['intersections'].set_ylabel('Substitutions in subset')
ax_dict['shading'].set_ylabel('Total per sample')
ax_dict['shading'].set_xlabel('Subset represented')

# %%
#Substitution mean Frequency by codon, destination

by = df.groupby(by=['Codon','Destination'])['logRatio'].mean().reset_index().pivot_table(columns='Destination',index='Codon',values='logRatio')
byann = df.groupby(by=['Codon','Destination'])['logRatio'].count().reset_index().pivot_table(columns='Destination',index='Codon',values='logRatio')
fig, ax = plt.subplots(figsize=(8,16))
ax = sns.heatmap(data = by,cmap='plasma_r',cbar_kws={'label': 'log10 Substitution Frequency'}, annot=byann, square=True, fmt='.0f', ax=ax)
ax.set_facecolor((.2,.2,.2,.15))
plt.savefig(os.path.join(dist_path,'Frequency mean by Codon and Destination'))

# %%
