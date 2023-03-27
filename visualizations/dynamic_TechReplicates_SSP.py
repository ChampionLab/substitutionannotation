#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:08:53 2022
A script to create many interactive plotly charts
Currently compatible with Fragger Output
Follows FindSubs_Fragpipe_PSM.py
@author: taylo
"""

#%%
#Try to make a plotly plot
import pandas as pd
import numpy as np
import plotly
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, html, dcc, Input, Output
import os
from substitutionannotation.resources import assets as p
import re

import upsetplot
import matplotlib.pyplot as plt
  
#%%
from io import BytesIO
import base64
def fig_to_uri(in_fig, close_all=True, **save_args):
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)
#%%
#Parse data for plotting

dfSubs = pd.read_csv(os.path.join(p.outputdir,'SSP Quant Tech.csv'))
dfAnnotate = pd.read_csv(os.path.join(p.outputdir,'SSP PSM.csv'))
#Subs with intensity, base peptide intensity, and protein intensity in the same sample,
# have already been processed by FindSubs_Fragpipe_PSM.py

missingoffset = 0.9

#Get unique Sample and replicate names
samples = dfSubs['Sample'].unique()
replicates = dfSubs['Replicate'].unique()

bySample = dfSubs.groupby(by=['Sample'])

#Define dictionary of Sample/Tech Rep combos
sampleReplicateDict = {}
for sample in samples:
    sampleReplicateDict[sample] = bySample.get_group(sample)['Replicate'].drop_duplicates().to_list()

     

#get FASTA and protein sequence info
fastadir = r"C:\Users\taylo\Documents\Notre Dame\Lab\FASTA\2022-03-26-decoys-contam-UP_2022_03_25_EcoliK12.fasta.fas"
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

dfAnnotate['Substitution Position'] =  dfAnnotate['Modified Peptide'].apply(
    lambda x:re.search('.\[',x).start()) + dfAnnotate['Protein Start']
#%%
#Start Dash

plotly.io.renderers.default='browser'
app = Dash(__name__,suppress_callback_exceptions=True)

#%%
#Plotly app stuff for interactive graphs

app.layout = html.Div(children=[
    dcc.RadioItems(options=[{'label':'Normalized Values','value':'logNormalized'},
                            {'label':'Ratio','value':'logRatio'},
                            {'label':'Intensity','value':'logIntensity'}],
                   id='inSignal', value = 'logNormalized'),
    html.Label('Select sample'),
    dcc.Dropdown(id='sampleSelected',options=samples,
                    value=samples[0]),
    html.Div(
        [html.Img(id = 'upsetPlot', src = '')],
        id='upsetPlot_div'),
    html.Label('Number of histogram bins'),
    dcc.Slider(id='nbins2',value=200, min=20, max=500),
    dcc.Graph(id='residualHistogram'),
    html.Br(),
    dcc.Dropdown(id='techRep1'),
    html.Br(),
    dcc.Dropdown(id='techRep2'),
    dcc.Graph(id='scatterTechReps')
],style={'padding':5,'flex':1})

#%%
@app.callback(
    Output('techRep1', 'options'),
    Output('techRep2','options'),
    Input('sampleSelected', 'value'))
def set_techRep_options(sample):
    return [{'label': i, 'value': i} for i in sampleReplicateDict[sample]],[{'label': i, 'value': i} for i in sampleReplicateDict[sample]]

#%%
@app.callback(
    Output('techRep1', 'value'),
    Input('techRep1', 'options')
    )
def set_techRep1(available_options):
    return available_options[0]['value']

#%%
@app.callback(
    Output('techRep2', 'value'),
    Input('techRep2', 'options')
    )
def set_techRep2(available_options):
    return available_options[1]['value']


#%%
@app.callback(
    Output('upsetPlot','src'),
    Output('residualHistogram','figure'),
    Input('sampleSelected','value'),
    Input('nbins2','value'),
    Input('inSignal','value'),
    )
def update_upsetPlot_and_redisualHistogramTechReps(sampleSelected,nbins2,inSignal):
     
    #UpSet plot things
    #Format input data correctly
    dict_upsetpep = {}
    dff = bySample.get_group(sampleSelected)
    peptide = 'Modified Peptide'
    dict_upsetpep = dff.groupby(by='Replicate')[peptide].apply(lambda x:
                            list(x.drop_duplicates())).to_dict()
            
    upsetpep = upsetplot.from_contents(dict_upsetpep)
    
    #Plot Upset plot
    figUpsetMPL,ax = plt.subplots()
    upsetplot.UpSet(upsetpep, subset_size='count',
                        show_percentages=True,
                        sort_categories_by=None).plot(fig=figUpsetMPL)
    plt.suptitle('Upset Plot of '+sampleSelected+' Unique Peptide Sequences')
    figUpset = fig_to_uri(figUpsetMPL)
    

        
    #Residual Histogram Parse
    dff = dff.groupby(by=['Replicate',peptide])[inSignal].max()
    dff = dff.reset_index().pivot(index=peptide,columns='Replicate',values=inSignal)
    #Remove missing values, impugned values
    dff.replace(0,np.nan,inplace=True)
 
    cv = pd.DataFrame()
    cv['cv'] = dff.apply(p.log10CV,axis=1)
    cv['n'] = dff.count(axis=1)
    #Remove peptides only observed in 1 replicate
    cv = cv[cv['n']>1]
    #Shove outliers to a specific extreme value
    cv.loc[cv['cv']>200,'cv'] = 200
    #Residual Histogram Plot
    figResidual = px.histogram(cv,x = 'cv', color = 'n', barmode = 'overlay',
                               nbins=nbins2,labels={'x':'%CV','y':'# of Peptides'})

    return figUpset, figResidual

#%%
@app.callback(
    Output('scatterTechReps','figure'),
    Input('sampleSelected','value'),
    Input('techRep1','value'),
    Input('techRep2','value'),
    Input('inSignal','value'),
)
def update_scatterTechReps(sampleSelected,techRep1,techRep2,inSignal):
    dff = bySample.get_group(sampleSelected).groupby(by=['Replicate'])
    peptide = 'Modified Peptide'


    xdata = dff.get_group(techRep1).groupby(by=peptide)[['Protein',inSignal]].max()
    xdata.columns = ['Protein',str(techRep1) +" "+inSignal]
    ydata = dff.get_group(techRep2).groupby(by=peptide)[['Protein',inSignal]].max()
    ydata.columns = ['Protein',str(techRep2) +" "+inSignal]
    dfplot = pd.merge(xdata,ydata,left_index=True,right_index=True,how='outer') 
    dfplot['Protein'] = dfplot['Protein_x'].fillna(dfplot['Protein_y'])
    dfplot.drop(columns=['Protein_x','Protein_y'],inplace=True)
    dfplot.reset_index(inplace=True)
    #Handle missing values

   #Mask for substitutions missing an intensity in one or more samples
    missingInX = dfplot[str(techRep1) +" "+inSignal].isna()
    missingInY = dfplot[str(techRep2) +" "+inSignal].isna()    
    notMissing = ~(missingInX | missingInY)
    dfplot.loc[missingInX,str(techRep1)+" "+inSignal] = dfplot[str(techRep1)+" "+inSignal].min()*missingoffset
    dfplot.loc[missingInY,str(techRep2)+" "+inSignal] = dfplot[str(techRep2)+" "+inSignal].min()*missingoffset
    
    #Compute minimum and maximum values to get range for plot
    xmin = dfplot[str(techRep1) +" "+inSignal].min()
    xmax = dfplot[str(techRep1) +" "+inSignal].max()
    ymin = dfplot[str(techRep2) +" "+inSignal].min()
    ymax = dfplot[str(techRep2) +" "+inSignal].max()
    
    #Scatter plot to compare intensities, frequencies at peptide level
    fig = plotly.subplots.make_subplots(rows=2, cols=2,
                    row_heights=[0.2,0.80],
                    column_widths=[0.80,0.2],
                    vertical_spacing = 0.005,
                    horizontal_spacing=0.005,
                    shared_yaxes=True,
                    shared_xaxes=True,
                    )
    
    #Main Figure
    fig.add_trace(go.Scatter(x=dfplot.loc[notMissing,str(techRep1) +" "+inSignal], 
                             y=dfplot.loc[notMissing,str(techRep2) +" "+inSignal],
                    hoverinfo='text',
                    mode='markers',
                    marker={'size':3},
                    marker_color='blue',
                    text=(dfplot[peptide]+'<br>'+dfplot['Protein']),
                    name = 'Shared Peptides'
                    ),
                    row=2,col=1)
    #y=x line
    fig.add_trace(go.Scatter(x=[-12,10],y=[-12,10],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1),
                            name='y=x'),
                            row=2,col=1)
    #20% above y=x line
    fig.add_trace(go.Scatter(x=[-12.079,10.079],y=[-12.079,10.079],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1),
                            name='y=x'),
                            row=2,col=1)
    #20% below y=x line
    fig.add_trace(go.Scatter(x=[-11.903,9.903],y=[-11.903,9.903],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1),
                            name='y=x'),
                            row=2,col=1)
    
    #Marginal plots
    fig.add_trace(go.Violin(x=dfplot.loc[missingInY,str(techRep1)+" "+inSignal],
                            side='negative',
                            box_visible=False,
                            showlegend=True,
                            marker={'opacity':0},
                            line_color='mediumpurple',
                            name='Unique to '+ str(techRep1),
                            y0=0,
                            scalegroup='right'),
                            row=1,col=1)
    fig.add_trace(go.Violin(y=dfplot.loc[missingInX,str(techRep2)+" "+inSignal],
                            side='negative',
                            box_visible=False,
                            showlegend=True,
                            marker={'opacity':0},
                            line_color='orange',
                            name='Unique to '+str(techRep2),
                            x0=0,
                            scalegroup='top'),
                            row=2,col=2)
    fig.add_trace(go.Violin(y=dfplot.loc[notMissing,str(techRep2)+" "+inSignal],
                            side='positive',
                            box_visible=False,
                            showlegend=False,
                            marker={'opacity':0},
                            line_color='blue',
                            name='Not Missing Violin',
                            x0=0,
                            scalegroup='top'),
                            row=2,col=2)
    fig.add_trace(go.Violin(x=dfplot.loc[notMissing,str(techRep1)+" "+inSignal],
                            side='positive',
                            box_visible=False,
                            showlegend=False,
                            marker={'opacity':0},
                            name='Not Missing Violin',
                            line_color='blue',
                            y0=0,
                            scalegroup='right'),
                            row=1,col=1)
    #Labels
    fig.update_layout(

        legend=dict(
            orientation='h',
            yanchor="bottom",
            y=1.01,
            xanchor="right",
            x=1,
            font=dict(size=12),
            bgcolor='rgba(0,0,0,0)',
            itemsizing= 'constant'
            ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        width=550,
        height=550)
    fig.update_xaxes(showgrid=True, gridwidth=0.125, gridcolor='rgba(200,200,200,66)',
                     zerolinewidth = 0.125, zerolinecolor='rgba(200,200,200,66)')
    fig.update_yaxes(showgrid=True, gridwidth=0.125, gridcolor='rgba(200,200,200,66)',
                     zerolinewidth= 0.125, zerolinecolor='rgba(200,200,200,66)')
    fig.update_xaxes(title_text=str(techRep1), range=[xmin,xmax],
                     tickmode = 'array',
                     tickvals = [x for x in range(-8,9,1)],
                     ticktext = [str(10**x) for x in range(-8,9,1)],
                     row=2, col=1)
    fig.update_yaxes(title_text=str(techRep2), range=[ymin,ymax],
                     tickmode = 'array',
                     tickvals = [x for x in range(-8,9,1)],
                     ticktext = [str(10**x) for x in range(-8,9,1)],
                     row=2, col=1)
    fig.update_xaxes(title_text='Density',showticklabels=False,row=2,col=2)
    fig.update_yaxes(title_text='Density',showticklabels=False,row=1,col=1)
    
    return fig
#%%
if __name__ == '__main__':
    app.run_server(debug=True)