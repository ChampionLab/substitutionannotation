# -*- coding: utf-8 -*-
"""
Created on  Oct 25 2023
A script to create many static pairwise comparison charts
Must follow FindSubs_xxx.py
@author: taylo
"""

#%%
#Import packages
import pandas as pd
import numpy as np
import seaborn as sns
import os
from substitutionannotation.resources import assets as sA
import re
import scipy.stats as st
import random
from matplotlib import pyplot as plt


#%%
#Define functions
def concatenate_lists(series):
    concatenated_list = []
    for lst in series:
        concatenated_list.extend(lst)
    return concatenated_list

def getCV(series):
    try:
        result = np.abs(np.std(series)/np.mean(series))
    except ZeroDivisionError:
        result = np.nan
    return result


#%%
#Parse data for plotting
missingoffset = 1.10

#Use filtered data if present. Otherwise, use generic output.
if os.path.exists(os.path.join(sA.outputdir,'Results Removed by Filter')):
    dfSubs = pd.read_csv(os.path.join(sA.outputdir,'Filtered SSP Quant.csv'))
    dfReplicates = pd.read_csv(os.path.join(sA.outputdir,'Filtered SSP Quant Replicates.csv'))   
else:
    dfSubs = pd.read_csv(os.path.join(sA.outputdir,'SSP Quant.csv'))
    dfReplicates = pd.read_csv(os.path.join(sA.outputdir,'SSP Quant Replicates.csv'))

dfAnnotate = pd.read_csv(os.path.join(sA.outputdir,'SSP PSM.csv'))
dfBaseBySample = pd.read_csv(os.path.join(sA.outputdir,'Base Peptides.csv'))

#Subs with intensity, with base peptide intensity and protein intensity
# in the same sample have already been processed by FindSubs_xxx.py

#Preparse some data
dfSubs = dfSubs.merge(dfAnnotate[['Modified Peptide','PTMs','Substitution Position']].drop_duplicates(),on='Modified Peptide',how='left')
dfReplicates = dfReplicates.merge(dfAnnotate[['Modified Peptide','PTMs']].drop_duplicates(),on='Modified Peptide',how='left')
samples = dfSubs['Sample'].drop_duplicates().reset_index(drop=True)
bySample = dfSubs.groupby(by=['Sample'])


dfReplicates['Ratio'] = dfReplicates['Intensity']/dfReplicates['Base Intensity']
dfReplicates['logRatio'] = np.log10(dfReplicates['Ratio'].replace(0,np.nan))



repBySample = dfReplicates.groupby(by=['Sample'])
bySamplePep = dfSubs.groupby(by=['Sample','Modified Peptide'])
dfBaseBySample['logIntensity'] = np.log10(dfBaseBySample['Base Intensity'].replace(0,np.nan))
#Determine biological relevance threshold
CVs = dfReplicates.groupby(by=['Sample','Modified Peptide'])['Intensity'].apply(getCV)
try:
    quantiles = pd.read_csv(os.path.join(sA.outputdir,'CV95 Quantiles.csv'))
    ratioCV = quantiles['Ratio CV'].item()
except:
    ratioCV = 1

#get FASTA and protein sequence info
# Get the directory path of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Get the parent directory of the script directory
module_dir = os.path.dirname(script_dir)
fastadir = os.path.join(module_dir,'resources','2022-03-26-decoys-contam-UP_2022_03_25_EcoliK12.fasta.fas')
file = open(fastadir,mode='r')
dfproteins = pd.DataFrame(columns=['Header','Sequence'])
i=1
accession = 'Failed'
while True:
    i +=1
    line = file.readline()
    if line == '':
        break
    if i%2 == 0:
        accession = re.search(r'\|(.*)\|',line).group(1)
        dfproteins.loc[accession,'Header'] = line
    if i%2 == 1:
        dfproteins.loc[accession,'Sequence'] = line

eftu = dfSubs.loc[dfSubs['Protein'].str.contains('EFTU'),'Protein'].max()



#define pairs for comparison
pairs = []
# Iterate through the list
for i in range(len(samples)):
    for j in range(i + 1, len(samples)):  # Avoid pairing an element with itself and duplicates
        pair = (samples[i], samples[j])
        pairs.append(pair)

#Check/make directory for outputs
output = os.path.join(sA.outputdir,'Static Pairwise Visualizations')
if not os.path.exists(output):
    os.mkdir(output)
for pair in pairs:
    if not os.path.exists(os.path.join(output,f'{pair}')):
        os.mkdir(os.path.join(output,f'{pair}'))

signaltype = 'logRatio'
signaltypeVolcano = 'Ratio'


#%%
#Make scatter plots
for pair in pairs:
    xdata = bySample.get_group(pair[0]).groupby(by='Modified Peptide')[['Protein',signaltype]].max().dropna(subset=[signaltype])
    xSignal = str(pair[0]) + ' ' + signaltype
    ySignal = str(pair[1]) + ' '+ signaltype
    xdata.columns = ['Protein',xSignal]
    ydata = bySample.get_group(pair[1]).groupby(by='Modified Peptide')[['Protein',signaltype]].max().dropna(subset=[signaltype])
    ydata.columns = ['Protein',ySignal]
    dfplot = pd.merge(xdata,ydata,left_index=True,right_index=True,how='outer') 
    dfplot['Protein'] = dfplot['Protein_x'].fillna(dfplot['Protein_y'])
    dfplot.drop(columns=['Protein_x','Protein_y'],inplace=True)
    dfplot.reset_index(inplace=True)
    
   #Mask for substitutions missing an intensity in one or more samples
    missingInX = dfplot[xSignal].isna()
    missingInY = dfplot[ySignal].isna()    

    notMissing = ~(missingInX | missingInY)
    dfplot.loc[missingInX,xSignal] = dfplot[xSignal].min()*missingoffset
    dfplot.loc[missingInY,ySignal] = dfplot[ySignal].min()*missingoffset

    #Min/max axis values
    xmin = dfplot.loc[notMissing,xSignal].min()
    ymin = dfplot.loc[notMissing,ySignal].min()
    xmax = dfplot.loc[notMissing,xSignal].max()
    ymax = dfplot.loc[notMissing,ySignal].max()
    squaremin = min([xmin,ymin])
    squaremax = max([xmax,ymax])

    #Compute p value of shared data, that x < y
    (t,pval) = st.ttest_rel(dfplot.loc[notMissing,xSignal],
                         dfplot.loc[notMissing,ySignal],
                         alternative ='less')
    
    #Scatter plot to compare intensities, frequencies at peptide level
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,
                    height_ratios=[0.2,0.80],
                    width_ratios=[0.80,0.2],
                    sharex='col',
                    sharey='row',
                    )
    plt.subplots_adjust(wspace=0.05, hspace=0.05)  
    #Main Figure
    ax3.scatter(x=dfplot.loc[notMissing,xSignal], 
                             y=dfplot.loc[notMissing,ySignal],
                    )
    # Y=X line
    ax3.plot([squaremin,squaremax],[squaremin,squaremax],
             linestyle='--',color='k')
    
    #Marginal plots
    sns.kdeplot(y=dfplot.loc[notMissing,ySignal],ax=ax4,label='Shared',
                common_norm=False, cut=0, color='blue' )
    sns.kdeplot(x=dfplot.loc[notMissing,xSignal],ax=ax1,
                common_norm=False, cut=0, color='blue' )
    sns.kdeplot(x=dfplot.loc[missingInY,xSignal],ax=ax1,label='Unique to '+str(pair[0]),
                common_norm=False, cut=0, color='purple' )
    sns.kdeplot(y=dfplot.loc[missingInX,ySignal],ax=ax4,label='Unique to '+str(pair[1]),
                common_norm=False, cut=0, color='orange' )
    
    #Labels
    ax3.set_ylabel(f'{pair[1]} substitution frequency')
    ax3.set_xlabel(f'{pair[0]} substitution frequency')
    ax2.set_axis_off()

    #Legend
    handles, labels = [], []
    for ax in (ax1, ax3, ax4):
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

    legend = fig.legend(handles, labels, loc=(0.75,0.775))
    fig.savefig(os.path.join(output,f'{pair}','Scatterplot'))

#%% Volcano plot
for pair in pairs:
    volcanoNumerator = pair[0]
    volcanoDenominator = pair[1]
    dff = pd.concat([repBySample.get_group(name) for name in [pair[0],pair[1]]])
    bySampleMetric = dff.groupby(by=['Sample','Modified Peptide'])[signaltypeVolcano].apply(list)
    bySampleMetric = bySampleMetric.reset_index().pivot_table(columns='Sample',index='Modified Peptide',values=signaltypeVolcano,
                                                     aggfunc=lambda x:x)
    #Find 0/Na in numerator or denominator for offsetting in plot
    maskNaNumerator = bySampleMetric[volcanoNumerator].isna()
    maskNaDenominator = bySampleMetric[volcanoDenominator].isna()
    maskBoth = ~(maskNaDenominator | maskNaNumerator)
    #Perform Mann-Whitney test, as we do not assume normality in these data
    bySampleMetric['P'] = bySampleMetric.loc[maskBoth].apply(lambda x:
                        st.ttest_ind(*x)[1],axis=1)
    
    #Calculate FC based on medians, adding filters for 0's.
    bySampleMetric['Fraction'] = (bySampleMetric[volcanoNumerator].apply(np.median)/
                                    bySampleMetric[volcanoDenominator].apply(np.median))
    maskNaNumerator = (maskNaNumerator | (bySampleMetric['Fraction'] == 0))
    maskNaDenominator = (maskNaDenominator | (bySampleMetric['Fraction'] == np.inf))
    bySampleMetric.loc[maskNaNumerator,'Fraction'] = np.nan
    bySampleMetric.loc[maskNaDenominator,'Fraction'] = np.nan
    bySampleMetric['FC'] = np.log2(bySampleMetric['Fraction'])
    
    #Perform Benjamini-Hochberg p-value correction 
    alpha = 0.05
    ntests=bySampleMetric['P'].count()    
    bySampleMetric['P Rank'] = bySampleMetric['P'].rank()
    bySampleMetric['Crit Value'] = (bySampleMetric['P Rank']/ntests)*alpha
    maskBen = bySampleMetric['P'] <= bySampleMetric['Crit Value']
    benpv= bySampleMetric[maskBen]['P'].max()
    if np.isnan(benpv):
        benpv = 0.05
   
    #Create mask for significant values
    maskSignificant = bySampleMetric['P'] <= benpv

    #Calculate biologically relevant cutoff
    if ratioCV >= 1: #If ratioCV >=100%, FC could approach 0 so it doesn't make sense to use CV based filter
        bioRel  = 1   #Use default 2-fold change cutoff
    else:
        #Get the distribution of fold change error
        dfbioRel = bySample.get_group(volcanoNumerator)[['Modified Peptide','Ratio CV']].merge(
            bySample.get_group(volcanoDenominator)[['Modified Peptide','Ratio CV']],how='outer',on='Modified Peptide')
        dfbioRel['Propigated CV'] = np.sqrt(dfbioRel['Ratio CV_x']**2+dfbioRel['Ratio CV_y']**2)
        propRatioCV = np.quantile(dfbioRel['Propigated CV'].dropna(),q=.95)
        if propRatioCV >= 1: #If propRatioCV >=100%, FC could approach 0 so it doesn't make sense to use CV based filter
            bioRel = 1
        else:
            bioRel = np.log2(1+propRatioCV)

    #Mask significant values AND biologically relevent values
    maskBioRel = (bySampleMetric['FC'] > bioRel) | (bySampleMetric['FC'] < -1*bioRel)
    maskSignificant = (maskSignificant & maskBioRel)


    bySampleMetric['logP'] = np.log10(bySampleMetric['P'])*-1
    #Get upper/lower bounds of data for plotting unique subset of data
    if (bySampleMetric['logP'] == np.inf).any():
        sigmax = bySampleMetric['logP'].replace(np.inf,np.nan).max()
        bySampleMetric.loc[bySampleMetric['logP'] == np.inf,'logP'] = sigmax*1.1
    if benpv == 0:
        benpv = 10**(bySampleMetric['logP'].max()*-1)
    sigmin = bySampleMetric['logP'].min()
    sigmax = bySampleMetric['logP'].max()
    fcmin = bySampleMetric['FC'].min()*missingoffset   
    fcmax = bySampleMetric['FC'].max()*missingoffset

    fig,ax = plt.subplots()
    #Add significantly enriched values
    ax.scatter(x=bySampleMetric.loc[maskSignificant & maskBoth,'FC'],
                                y=bySampleMetric.loc[maskSignificant & maskBoth,'logP'],
                            color='blue', label='Significant')
    #Add insignificant values
    ax.scatter(x=bySampleMetric.loc[~maskSignificant & maskBoth,'FC'],
                                y=bySampleMetric.loc[~maskSignificant & maskBoth,'logP'],
                                color = 'grey')

    #Add values only in one sample
    #random x is to add jitter
    if maskNaNumerator.any():
        ax.scatter(x=[fcmin*(1+random.randrange(1,10)/100) for x in range(len(bySampleMetric.loc[maskNaNumerator]))],
                                y=np.linspace(sigmin,sigmax,len(bySampleMetric.loc[maskNaNumerator])),
                                color='purple', label='Unique to ' + str(volcanoDenominator))
    if maskNaDenominator.any():
        ax.scatter(x=[fcmax*(1+random.randrange(1,10)/100) for x in range(len(bySampleMetric.loc[maskNaDenominator]))],
                                y=np.linspace(sigmin,sigmax,len(bySampleMetric.loc[maskNaDenominator])),
                                color='orange', label='Unique to '+str(volcanoNumerator))
    #Add dotted lines
    #Significance Threshold
    ax.plot([fcmin*1.2,fcmax*1.2],[-1*np.log10(benpv),-1*np.log10(benpv)],
                                linestyle='--', color='k')
    #FC lower threshold
    ax.plot([-bioRel,-bioRel],[0,sigmax*1.2],linestyle='--', color='k')
    #FC upper threshold
    ax.plot([bioRel,bioRel],[0,sigmax*1.2],linestyle='--', color='k')
    #Missing data upper fence
    ax.plot([fcmax*0.98,fcmax*0.98],[0,sigmax*1.2],linestyle=':', color='k')
    #Missing data lower fence
    ax.plot([fcmin*0.98,fcmin*0.98],[0,sigmax*1.2],linestyle=':', color='k')

    #Labels
    ax.set_xlabel(f'log2({volcanoNumerator}/{volcanoDenominator})')
    ax.set_yticks(ticks = [x for x in range(0,round(sigmax)+1)],
                labels = [str(10**-x) for x in range(0,round(sigmax)+1)])
    ax.set_ylabel('-10log(P)')
    ax.spines['top'].set_visible(False)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3)
    plt.savefig(os.path.join(output,f'{pair}','Volcano'))
   

   #Output table information
    summary = pd.Series(index = ['# Peptides Shared',' # Peptides Significant','# of Peptides Up','# of Peptides Down'
                                      ,f'# Unique to {volcanoNumerator}',f'# Unique to {volcanoDenominator}'])
    
    maskUp = maskSignificant & maskBoth & (bySampleMetric['FC'] > 0)
    maskDown = maskSignificant & maskBoth & (bySampleMetric['FC'] < 0)
    summary['# Peptides Shared'] = len(bySampleMetric[maskBoth])
    summary['# Peptides Significant'] = len(bySampleMetric[maskBoth & maskSignificant])
    summary['# of Peptides Up'] = len(bySampleMetric[maskUp])
    summary['# of Peptides Down'] = len(bySampleMetric[maskDown])
    summary[f'# Unique to {volcanoNumerator}'] = len(bySampleMetric[maskNaDenominator])
    summary[f'# Unique to {volcanoDenominator}'] = len(bySampleMetric[maskNaNumerator])
    summary.to_csv(os.path.join(output,f'{pair}',f'{pair} Summary.csv'))

# %%
