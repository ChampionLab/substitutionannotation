# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:08:53 2022
A script to create many interactive plotly charts
Must follow FindSubs_xxx.py
@author: taylo
"""

#%%
#Import packages
import pandas as pd
import numpy as np
import plotly
from plotly import subplots as psubplots
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, html, dcc, Input, Output
import os
from substitutionannotation.resources import assets as sA

import re
import scipy.stats as st
import random



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

dfSubs = pd.read_csv(os.path.join(sA.outputdir,'SSP Quant.csv'))
dfReplicates = pd.read_csv(os.path.join(sA.outputdir,'SSP Quant Tech.csv'))
dfAnnotate = pd.read_csv(os.path.join(sA.outputdir,'SSP PSM.csv'))
#Subs with intensity, with base peptide intensity and protein intensity
# in the same sample have already been processed by FindSubs_xxx.py

#Preparse some data
dfSubs = dfSubs.merge(dfAnnotate[['Modified Peptide','PTMs','Substitution Position']].drop_duplicates(),on='Modified Peptide',how='left')
dfReplicates = dfReplicates.merge(dfAnnotate[['Modified Peptide','PTMs']].drop_duplicates(),on='Modified Peptide',how='left')
samples = dfSubs['Sample'].drop_duplicates().reset_index(drop=True)
bySample = dfSubs.groupby(by=['Sample'])
repBySample = dfReplicates.groupby(by=['Sample'])
bySamplePep = dfSubs.groupby(by=['Sample','Modified Peptide'])
quantiles = pd.read_csv(os.path.join(sA.outputdir,'CV95 Quantiles.csv'))
#Determine biological relevance threshold
CVs = dfReplicates.groupby(by=['Sample','Modified Peptide'])['Intensity'].apply(getCV)
bioRel = quantiles['logRatio CV'].item()

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

#Check/make directory for outputs
output = os.path.join(sA.outputdir,'dynamicVisualizations')
if not os.path.exists(output):
    os.mkdir(output)

#%%
#Start Dash

plotly.io.renderers.default='browser'
app = Dash(__name__,suppress_callback_exceptions=True)

#%%
#Plotly app stuff for interactive graphs

app.layout = html.Div([
    dcc.Dropdown(options = ['logRatio','logNormalized'], value = 'logRatio', id='signaltype'),
    dcc.Tabs(id="tabs", value='tab-1', children=[
        dcc.Tab(label='Ratio by Type', value='tab-1'),
        dcc.Tab(label='Ratio Scatterplot', value='tab-2'),
        dcc.Tab(label='Substitution Frequency by Protein', value='tab-3'),
        dcc.Tab(label='Ratio Volcano Plot',value='tab-4')
    ]),
    html.Div(id='tabs-content')
])

@app.callback(Output('tabs-content', 'children'),
              Input('tabs', 'value'))
def render_content(tab):
    if tab == 'tab-1':
        return html.Div([
            html.Div(children=[
                html.Br(),
                html.Label('Substitution(Regex) Filter 1: '),
                dcc.Input(id='subfilter',value='Gly', type='text'),
                html.Br(),
                dcc.Graph(
                    id='filtered_freq_violin'),
                    html.Br(),
                dcc.Graph(id='typehex',
                          style={
                'width': 'auto',
                'height': 'auto',
                'max-width': '90vw',  # Set the maximum width of the plot area
                'max-height': '90vh',  # Set the maximum height of the plot area
                'padding': '10px',  # Add padding around the plot area
                'box-sizing': 'border-box',  # Include padding in the dimensions
            })
            ], style={'padding': 5, 'flex': 1})
        ])
    elif tab == 'tab-2':
        return html.Div([
            html.Div(children=[
                html.H3('Substitution Ratio by Sample'
                        ),
                dcc.Dropdown(options=samples,value=samples[0],id='freq_scatter_x'),
                dcc.Dropdown(options=samples,value=samples[1],id='freq_scatter_y'),
                dcc.Graph(
                    id='freq_scatter'),
                ], style={'padding': 5, 'flex': 1}),
            html.Br(),
            html.Div(id='ttestOutput')
        ], style={'display': 'flex', 'flex-direction': 'row'}
            )
    elif tab== 'tab-3':
        return html.Div([
            html.Div(children=[
                html.H3('Substitution Frequency by Protein'),
                html.Br(),
                html.Label('Samples to compare'),
                dcc.Dropdown(options=samples,value=samples[0],id='freq_prot_top'),
                dcc.Dropdown(options=samples,value=samples[1],id='freq_prot_bottom'),
                html.Br(),
                dcc.Dropdown(options=dfSubs['Protein'].unique(),id='dropProtein',
                             value=eftu),
                dcc.Graph(
                    id='notManhattan')
                ],style={'padding':5,'flex'
                         :1}),
                         ])
    elif tab== 'tab-4':
        return html.Div([
            html.Div(children=[
                html.H3('Volcano Plot'),
                html.Br(),
                dcc.Dropdown(id='volcanoType',options=['Modified Peptide','PTMs'],
                             value='PTMs'),
                html.Br(),
                dcc.Dropdown(id='volcanoNumerator',options=samples,value=samples[0]),
                html.Br(),
                dcc.Dropdown(id='volcanoDenominator',options=samples,value=samples[1]),
                dcc.Graph(
                    id='volcano')
                ],style={'padding':5,'flex'
                         :1}
                )
            ])
                      
    



#%%
@app.callback(
    Output(component_id='filtered_freq_violin', component_property='figure'),
    Output(component_id='typehex',component_property= 'figure'),
    Input(component_id='subfilter', component_property='value'),
    Input(component_id ='signaltype', component_property= 'value')
)
def update_freq_violin(subfilter,signaltype):
    #Parse proper substitution type
    dff = dfSubs[dfSubs['PTMs'].str.contains(subfilter)]
    dff = dff.groupby(by='Sample')

    # Scale the size of each violin plot based on the number of peptides
    # Create a list of scales for each sample based on the number of peptides
    scales = pd.Series(dtype='float64')
    for sample in samples:
        try:
            scales.loc[sample] = len(dff.get_group(sample))
        except:     #Handle error when a sample has no peptides that match the filter
            continue
    max_scale = scales.max()
    scales = scales/max_scale

    fig = go.Figure()
    for sample in samples:
        try:
            dfplot = dff.get_group(sample)[[signaltype,'PTMs']].dropna(subset=signaltype)
        except:     #Handle error when a sample has no peptides that match the filter
            continue
        fig.add_trace(
        go.Violin(x0=str(sample),
                    y=dfplot[signaltype],
                    box_visible=True,points='all',
                    hoverinfo='text',text=dfplot['PTMs'],
                    scalemode='count', scalegroup=1,
                  showlegend=False)
                )
    fig.update_yaxes(title_text=signaltype)
    fig.update_layout(violingap=0,showlegend=False,title=subfilter + ' Substitutions')

    """ Rectangle density heatmap
    fig2 = px.density_heatmap(dfSubs[dfSubs['PTMs'].str.contains(subfilter)], x='logIntensity',
                               y=signaltype, nbinsx=30, nbinsy=30, color_continuous_scale='Blues')
    """
    #Traditional scatter
    fig2 = px.scatter(data_frame= dfSubs[dfSubs['PTMs'].str.contains(subfilter)], x="logIntensity",
                       y=signaltype)
    fig2.update_traces(marker=dict(size=5))
    
    return fig,fig2


#%%
@app.callback(
    Output('freq_scatter','figure'),
    Output('ttestOutput','children'),
    Input('freq_scatter_x','value'),
    Input('freq_scatter_y','value'),
    Input(component_id ='signaltype', component_property= 'value')

)
def update_frequency_scatter(freq_scatter_x,freq_scatter_y,signaltype):
    
    xdata = bySample.get_group(freq_scatter_x).groupby(by='Modified Peptide')[['Protein',signaltype]].max().dropna(subset=[signaltype])
    xSignal = str(freq_scatter_x) + ' ' + signaltype
    ySignal = str(freq_scatter_y) + ' '+ signaltype
    xdata.columns = ['Protein',xSignal]
    ydata = bySample.get_group(freq_scatter_y).groupby(by='Modified Peptide')[['Protein',signaltype]].max().dropna(subset=[signaltype])
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
    squaremax = max([xmax,xmin])

    #Compute p value of shared data, that x < y
    (t,pval) = st.ttest_rel(dfplot.loc[notMissing,xSignal],
                         dfplot.loc[notMissing,ySignal],
                         alternative ='less')
    
    #Scatter plot to compare intensities, frequencies at peptide level
    fig = psubplots.make_subplots(rows=2, cols=2,
                    row_heights=[0.2,0.80],
                    column_widths=[0.80,0.2],
                    vertical_spacing = 0.005,
                    horizontal_spacing=0.005,
                    shared_yaxes=True,
                    shared_xaxes=True,
                    )
    
    #Main Figure
    fig.add_trace(go.Scatter(x=dfplot.loc[notMissing,xSignal], 
                             y=dfplot.loc[notMissing,ySignal],
                    hoverinfo='text',
                    mode='markers',
                    marker={'size':3},
                    marker_color='blue',
                    text=(dfplot['Modified Peptide']+'<br>'+dfplot['Protein']),
                    name = 'Shared Peptides'
                    ),
                    row=2,col=1)
    # Y=X line
    fig.add_trace(go.Scatter(x=[squaremin,squaremax],y=[squaremin,squaremax],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1),
                            name='y=x'),
                            row=2,col=1)
    
    #Marginal plots
    fig.add_trace(go.Violin(x=dfplot.loc[missingInY,xSignal],
                            side='negative',
                            box_visible=False,
                            showlegend=True,
                            marker={'opacity':0},
                            line_color='mediumpurple',
                            name='Unique to '+ str(freq_scatter_x),
                            y0=0,
                            scalegroup='right',scalemode='count'),
                            row=1,col=1)
    fig.add_trace(go.Violin(y=dfplot.loc[missingInX,ySignal],
                            side='negative',
                            box_visible=False,
                            showlegend=True,
                            marker={'opacity':0},
                            line_color='orange',
                            name='Unique to '+str(freq_scatter_y),
                            x0=0,
                            scalegroup='top',scalemode='count'),
                            row=2,col=2)
    fig.add_trace(go.Violin(y=dfplot.loc[notMissing,ySignal],
                            side='positive',
                            box_visible=False,
                            showlegend=False,
                            marker={'opacity':0},
                            line_color='blue',
                            name='Not Missing Violin',
                            x0=0,
                            scalegroup='top',scalemode='count'),
                            row=2,col=2)
    fig.add_trace(go.Violin(x=dfplot.loc[notMissing,xSignal],
                            side='positive',
                            box_visible=False,
                            showlegend=False,
                            marker={'opacity':0},
                            name='Not Missing Violin',
                            line_color='blue',
                            y0=0,
                            scalegroup='right',scalemode='count'),
                            row=1,col=1)
    #Labels
    fig.update_layout(
        #title={'text':"Comparing Substitution Frequency by Sample",
        #       'y':0.9,
        #       'x':0.5,
        #       'xanchor':'center','yanchor':'top'},
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
    fig.update_xaxes(title_text=str(freq_scatter_x), row=2, col=1)
    fig.update_yaxes(title_text=str(freq_scatter_y), row=2, col=1)
    fig.update_xaxes(title_text='Density',showticklabels=False,row=2,col=2)
    fig.update_yaxes(title_text='Density',showticklabels=False,row=1,col=1)
    
    textout = 'One-sided T-test that '+str(freq_scatter_y)+' > ' + str(freq_scatter_x) +' = ' +"{:.2e}".format(pval)
    return fig,textout
#%%

@app.callback(
    Output(component_id='notManhattan', component_property='figure'),
    Input(component_id='dropProtein', component_property='value'),
    Input(component_id='freq_prot_top', component_property='value'),
    Input(component_id='freq_prot_bottom', component_property='value'),
    Input(component_id='signaltype',component_property='value')
)
def updateNotManhattan(dropProtein,sampleTop,sampleBot,signaltype):
    dff = dfSubs[(dfSubs['Protein'] == dropProtein)&(dfSubs['Sample']==sampleTop)]
    accession = re.search('(?<=\|).*(?=\|)',dropProtein)[0]
    sequence = dfproteins.loc[accession,'Sequence']    
    dff = dff.groupby(by=['PTMs','Substitution Position'])[signaltype].max()
    dff = dff.reset_index()
    dff['Destination'] = dff['PTMs'].str.extract('->(...)')
    
    dff2 = dfSubs[(dfSubs['Protein'] == dropProtein)&(dfSubs['Sample']==sampleBot)]  
    dff2 = dff2.groupby(by=['PTMs','Substitution Position'])[signaltype].max()
    dff2 = dff2.reset_index()
    dff2[signaltype] = -1*dff2[signaltype]
    dff2['Destination'] = dff2['PTMs'].str.extract('->(...)')
    
    #Make a color mapping dictionary 
    color_mapping = {}
    num_new_colors = 0
    for dest in pd.concat([dff['Destination'],dff2['Destination']]).drop_duplicates():
            color_mapping[dest] = px.colors.qualitative.Dark24[num_new_colors % len(px.colors.qualitative.D3)]
            num_new_colors += 1

    #I don't know how to use the dictionary on the go.scatter, so get a list
    color_codes = []
    for dest in dff2['Destination']:
        color_codes.append(color_mapping[dest])

    fig = px.scatter(dff, x='Substitution Position', 
                                y=signaltype,
                                color='Destination',
                                hover_data={'PTMs':True,
                                            signaltype:False,
                                            'Substitution Position':False,
                                            'Destination':False},
                                color_discrete_map=color_mapping
                                )
    
    fig.add_trace(go.Scatter(x=dff2['Substitution Position'],
                          y=dff2[signaltype],
                          mode='markers',
                          marker=dict(color=color_codes,
                                      symbol='circle'),
                          hovertemplate='%{text}',
                          text=dff2['PTMs'])
                )
    
    fig.layout.update(showlegend=False)#{'title':'Substitution Destination'})
    fig.update_xaxes(title_text=dropProtein, range=[0,len(sequence)],
                     tickmode = 'array',
                     tickvals = [x for x in range(0,len(sequence))],
                     ticktext = [x for x in sequence],
                     tickangle=0
                     )
    fig.update_yaxes(title_text='Substitution Frequency',
                     tickmode = 'array',
                     tickvals = [x for x in range(-12,12,1)],
                     ticktext = [r'10<sup>-{}</sup>'.format(np.abs(x)) for x in range(-12,12,1)]
                     )
    return fig



#%%
@app.callback(
    Output('volcano','figure'),
    Input('volcanoType','value'),
    Input('volcanoNumerator','value'),
    Input('volcanoDenominator','value'),
    Input('signaltype','value')
    )
def updateVolcanoPlot(volcanoType,volcanoNumerator,volcanoDenominator,signaltype):
    #Parse relevant data
    dff = pd.concat([repBySample.get_group(name) for name in [volcanoNumerator,volcanoDenominator]])
    bySampleMetric = dff.groupby(by=['Sample',volcanoType])[signaltype].apply(list)
    bySampleMetric = bySampleMetric.reset_index().pivot_table(columns='Sample',index=volcanoType,values=signaltype,
                                                     aggfunc=lambda x:x)
    #Find 0/Na in numerator or denominator for offsetting in plot
    maskNaNumerator = bySampleMetric[volcanoNumerator].isna()
    maskNaDenominator = bySampleMetric[volcanoDenominator].isna()
    maskBoth = ~(maskNaDenominator | maskNaNumerator)
    #Perform Mann-Whitney test, as we do not assume normality in these data
    bySampleMetric['P'] = bySampleMetric.loc[maskBoth].apply(lambda x:
                        st.mannwhitneyu(*x)[1],axis=1)
    #Calculate FC based on medians, adding filters for 0's.
    bySampleMetric['Fraction'] = (bySampleMetric[volcanoNumerator].apply(np.median)/
                                    bySampleMetric[volcanoDenominator].apply(np.median))
    maskNaNumerator = (maskNaNumerator | (bySampleMetric['Fraction'] == 0))
    maskNaDenominator = (maskNaDenominator | (bySampleMetric['Fraction'] == np.inf))
    bySampleMetric.loc[maskNaNumerator,'Fraction'] = np.nan
    bySampleMetric.loc[maskNaDenominator,'Fraction'] = np.nan
    bySampleMetric['FC'] = np.log2(bySampleMetric['Fraction'])
    
    #Perform Benhamini-Hochberg p-value correction 
    alpha = 0.05
    ntests=bySampleMetric['P'].count()    
    bySampleMetric['P Rank'] = bySampleMetric['P'].rank()
    bySampleMetric['Crit Value'] = (bySampleMetric['P Rank']/ntests)*alpha
    maskBen = bySampleMetric['P'] <= bySampleMetric['Crit Value']
    benpv= bySampleMetric[maskBen]['P'].max()
    #Create mask for significant values
    maskSignificant = bySampleMetric['P'] <= benpv
    #Mask significant values AND biologically relevent values

    maskBioRel = (bySampleMetric['FC'] > bioRel) | (bySampleMetric['FC'] < -1*bioRel)
    maskSignificant = (maskSignificant & maskBioRel)


    bySampleMetric['logP'] = np.log10(bySampleMetric['P'])*-1
    #Get upper/lower bounds of data for plotting unique subset of data
    sigmin = bySampleMetric['logP'].min()
    sigmax = bySampleMetric['logP'].max()
    fcmin = bySampleMetric['FC'].min()*missingoffset   
    fcmax = bySampleMetric['FC'].max()*missingoffset

    fig = go.Figure()
    #Add significantly enriched values
    fig.add_trace(go.Scatter(x=bySampleMetric.loc[maskSignificant & maskBoth,'FC'],
                             y=bySampleMetric.loc[maskSignificant & maskBoth,'logP'],
                             mode='markers',marker_color='rgba(255,0,0)', marker = dict(size=6),
                            name='Significant',
                            hoverinfo='text',
                            text=bySampleMetric.loc[maskSignificant].index)
                  )
    #Add values only in one sample
    #random x is to add jitter
    if maskNaNumerator.any():
        fig.add_trace(go.Scatter(x=[fcmin*(1+random.randrange(1,10)/100) for x in range(len(bySampleMetric.loc[maskNaNumerator]))],
                             y=np.arange(sigmin,sigmax,sigmax/len(bySampleMetric.loc[maskNaNumerator])),
                             mode='markers',marker_color='mediumpurple', marker = dict(size=3),
                            name='Unique to ' + str(volcanoDenominator),
                            hoverinfo='text',
                            text=bySampleMetric.loc[maskNaNumerator].index)
                  )
    if maskNaDenominator.any():
        fig.add_trace(go.Scatter(x=[fcmax*(1+random.randrange(1,10)/100) for x in range(len(bySampleMetric.loc[maskNaDenominator]))],
                             y=np.arange(sigmin,sigmax,sigmax/len(bySampleMetric.loc[maskNaDenominator])),
                             mode='markers',marker_color='orange', marker = dict(size=3),
                            name='Unique to '+str(volcanoNumerator),
                            hoverinfo='text',
                            text=bySampleMetric.loc[maskNaDenominator].index)
                  )
    #Add insignificant values
    fig.add_trace(go.Scatter(x=bySampleMetric.loc[~maskSignificant & maskBoth,'FC'],
                             y=bySampleMetric.loc[~maskSignificant & maskBoth,'logP'],
                             mode='markers',marker_color='grey', marker = dict(size=4),
                            name='Insignificant',
                            hoverinfo='text',
                            text=bySampleMetric.index)
                  )
    #Add dotted lines
    #Significance Threshold
    fig.add_trace(go.Scatter(x=[fcmin*1.2,fcmax*1.2],y=[-1*np.log10(benpv),-1*np.log10(benpv)],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1,dash='dot'),
                            name='Significance line')
                  )
    #FC lower threshold
    if 1-bioRel > 0:
        fig.add_trace(go.Scatter(x=[-bioRel,-bioRel],y=[0,sigmax*1.2],
                                mode='lines',marker_color='rgba(130,130,130)',
                                showlegend=False,
                                line=dict(width=1,dash='dot'),
                                name='FC Neg line')
                    )
    #FC upper threshold
    fig.add_trace(go.Scatter(x=[bioRel,bioRel],y=[0,sigmax*1.2],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1,dash='dot'),
                            name='FC Pos line')
                  )
    #Missing data upper fence
    fig.add_trace(go.Scatter(x=[fcmax*0.98,fcmax*0.98],y=[0,sigmax*1.2],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1,dash='solid'),
                            name='Missing upper fence')
                  )
    #Missing data lower fence
    fig.add_trace(go.Scatter(x=[fcmin*0.98,fcmin*0.98],y=[0,sigmax*1.2],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1,dash='solid'),
                            name='Missing lower fence')
                  )
    #Text annotations
    fig.add_annotation(text="Higher "+'Ratio'+" in "+ str(volcanoNumerator),
                xref="paper", yref="paper",
                x=.85, y=-0.12,showarrow=False)
    fig.add_annotation(text="Higher "+'Ratio'+" in "+ str(volcanoDenominator),
                xref="paper", yref="paper",
                x=.15, y=-0.12,showarrow=False)
    
    fig.update_xaxes(title_text='log2(fold change)', range=[fcmin*1.2,fcmax*1.2],
                     
                     )
    fig.update_yaxes(title_text='P', range=[0,round(sigmax)+1],
                     tickmode = 'array',
                     tickvals = [x for x in range(0,round(sigmax)+1)],
                     ticktext = [str(10**-x) for x in range(0,round(sigmax)+1)]
                     )
    
    #Output significant, unique data
    filestr = f'VolcanoSignificant {signaltype} {volcanoNumerator}v{volcanoDenominator}.csv'
    outputcsv = os.path.join(sA.outputdir,'dynamicVisualizations',filestr)

    if os.path.exists(outputcsv):#No need to remake the output every time the app refreshes
        pass
    else:
        bySampleMetric.loc[maskSignificant & maskBoth].to_csv(outputcsv)
    return fig

#%%
if __name__ == '__main__':
    app.run(debug=True)