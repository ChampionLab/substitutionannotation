# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:08:53 2022
A script to create many interactive plotly charts
Follows MSFraggerVis4
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
import scipy.stats as st
import random
import upsetplot
import matplotlib.pyplot as plt
import seaborn as sns
#%%
def concatenate_lists(series):
    concatenated_list = []
    for lst in series:
        concatenated_list.extend(lst)
    return concatenated_list
#%%
from io import BytesIO
import base64
def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
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
missingoffset = 1.10

dfSubs = pd.read_csv(os.path.join(p.outputdir,'SSP Quant.csv'))
dfAnnotate = pd.read_csv(os.path.join(p.outputdir,'SSP PSM.csv'))
#Subs with intensity, with base peptide intensity and protein intensity
# in the same sample have already been processed by FindSubs_Pragpipe_PSM.py


dfSubs = dfSubs.merge(dfAnnotate[['Modified Peptide','PTMs']],on='Modified Peptide',how='left')
samples = dfSubs['Sample'].drop_duplicates().reset_index(drop=True)
bySample = dfSubs.groupby(by=['Sample'])
bySamplePep = dfSubs.groupby(by=['Sample','Modified Peptide'])

                        



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

#%%
#Start Dash

plotly.io.renderers.default='browser'
app = Dash(__name__,suppress_callback_exceptions=True)

#%%
#Plotly app stuff for interactive graphs

app.layout = html.Div([
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
                    id='filtered_freq_violin')
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
                html.Label('Multi-Select Dropdown'),
                dcc.Dropdown(options=samples,value=samples[0],id='freq_prot_top'),
                dcc.Dropdown(options=samples,value=samples[1],id='freq_prot_bottom'),
                html.Br(),
                dcc.Dropdown(options=dfSubs['Protein'].unique(),id='dropProtein',
                             value='P0CE47|EFTU1'),
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
    Input(component_id='subfilter', component_property='value'),

)
def update_freq_violin(subfilter):
    #Parse proper substitution type
    dff = dfSubs[dfSubs['PTMs'].str.contains(subfilter)]
    dff = dff.groupby(by='Sample')

    fig = go.Figure()
    for sample in samples:
        try:
            dfplot = dff.get_group(sample)[['logRatio','PTMs']].dropna(subset='logRatio')
        except:
            continue
        fig.add_trace(
        go.Violin(x0=str(sample),
                    y=dfplot['logRatio'],
                    box_visible=True,points='all',
                    hoverinfo='text',text=dfplot['PTMs'],
                  showlegend=False)
                )
    fig.update_yaxes(title_text='log10 Ratio')
    fig.update_layout(violingap=0,showlegend=False,title=subfilter + ' Substitutions')
    return fig


#%%
@app.callback(
    Output('freq_scatter','figure'),
    Output('ttestOutput','children'),
    Input('freq_scatter_x','value'),
    Input('freq_scatter_y','value')
)
def update_frequency_scatter(freq_scatter_x,freq_scatter_y):
    
    xdata = bySample.get_group(freq_scatter_x).groupby(by='Modified Peptide')[['Protein','logNormalized']].max()
    xdata.columns = ['Protein',freq_scatter_x +" logNormalized"]
    ydata = bySample.get_group(freq_scatter_y).groupby(by='Modified Peptide')[['Protein','logNormalized']].max()
    ydata.columns = ['Protein',freq_scatter_y +" logNormalized"]
    dfplot = pd.merge(xdata,ydata,left_index=True,right_index=True,how='outer') 
    dfplot['Protein'] = dfplot['Protein_x'].fillna(dfplot['Protein_y'])
    dfplot.drop(columns=['Protein_x','Protein_y'],inplace=True)
    dfplot.reset_index(inplace=True)
    
   #Mask for substitutions missing an intensity in one or more samples
    missingInX = dfplot[freq_scatter_x +" logNormalized"].isna()
    missingInY = dfplot[freq_scatter_y +" logNormalized"].isna()    

    notMissing = ~(missingInX | missingInY)
    dfplot.loc[missingInX,freq_scatter_x+" logNormalized"] = dfplot[freq_scatter_x+" logNormalized"].min()*missingoffset
    dfplot.loc[missingInY,freq_scatter_y+" logNormalized"] = dfplot[freq_scatter_y+" logNormalized"].min()*missingoffset
    #Compute p value of shared data, that x < y
    (t,pval) = st.ttest_rel(dfplot.loc[notMissing,freq_scatter_x +" logNormalized"],
                         dfplot.loc[notMissing,freq_scatter_y +" logNormalized"],
                         alternative ='less')
    
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
    fig.add_trace(go.Scatter(x=dfplot.loc[notMissing,freq_scatter_x +" logNormalized"], 
                             y=dfplot.loc[notMissing,freq_scatter_y +" logNormalized"],
                    hoverinfo='text',
                    mode='markers',
                    marker={'size':3},
                    marker_color='blue',
                    text=(dfplot['Modified Peptide']+'<br>'+dfplot['Protein']),
                    name = 'Shared Peptides'
                    ),
                    row=2,col=1)

    fig.add_trace(go.Scatter(x=[-15,-3],y=[-15,-3],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1),
                            name='y=x'),
                            row=2,col=1)
    
    #Marginal plots
    fig.add_trace(go.Violin(x=dfplot.loc[missingInY,freq_scatter_x+" logNormalized"],
                            side='negative',
                            box_visible=False,
                            showlegend=True,
                            marker={'opacity':0},
                            line_color='mediumpurple',
                            name='Unique to '+ freq_scatter_x,
                            y0=0,
                            scalegroup='right'),
                            row=1,col=1)
    fig.add_trace(go.Violin(y=dfplot.loc[missingInX,freq_scatter_y+" logNormalized"],
                            side='negative',
                            box_visible=False,
                            showlegend=True,
                            marker={'opacity':0},
                            line_color='orange',
                            name='Unique to '+freq_scatter_y,
                            x0=0,
                            scalegroup='top'),
                            row=2,col=2)
    fig.add_trace(go.Violin(y=dfplot.loc[notMissing,freq_scatter_y+" logNormalized"],
                            side='positive',
                            box_visible=False,
                            showlegend=False,
                            marker={'opacity':0},
                            line_color='blue',
                            name='Not Missing Violin',
                            x0=0,
                            scalegroup='top'),
                            row=2,col=2)
    fig.add_trace(go.Violin(x=dfplot.loc[notMissing,freq_scatter_x+" logNormalized"],
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
    fig.update_xaxes(title_text=freq_scatter_x, row=2, col=1)
    fig.update_yaxes(title_text=freq_scatter_y, row=2, col=1)
    fig.update_xaxes(title_text='Density',showticklabels=False,row=2,col=2)
    fig.update_yaxes(title_text='Density',showticklabels=False,row=1,col=1)
    
    textout = 'One-sided T-test that '+freq_scatter_y+' > ' + freq_scatter_x +' = ' +"{:.2e}".format(pval)
    return fig,textout
#%%

@app.callback(
    Output(component_id='notManhattan', component_property='figure'),
    Input(component_id='dropProtein', component_property='value'),
    Input(component_id='dropSample', component_property='value')
)
def updateNotManhattan(dropProtein,dropSample):
    dff = dfSubs[(dfSubs['Protein'].str.contains(dropProtein)&dfSubs['Sample'].isin(dropSample))]
    accession = re.match('.*(?=\|)',dropProtein)[0]
    sequence = dfproteins.loc[accession,'Sequence']    
    dff = dff.groupby(by=['PTMs','Substitution Position'])['logNormalized'].max()
    dff = dff.reset_index()
    dff['Destination'] = dff['PTMs'].str.extract('->(...)')
    
    fig = px.scatter(dff, x='Substitution Position', 
                                y="logNormalized",
                                color='Destination',
                                hover_data={'PTMs':True,
                                            'logNormalized':False,
                                            'Substitution Position':False,
                                            'Destination':False},
                                color_discrete_sequence=px.colors.qualitative.Dark24,
                                )
    
    fig.layout.update(showlegend=False)#{'title':'Substitution Destination'})
    fig.update_xaxes(title_text=dropProtein, range=[0,len(sequence)],
                     tickmode = 'array',
                     tickvals = [x for x in range(0,len(sequence))],
                     ticktext = [x for x in sequence],
                     tickangle=0
                     )
    fig.update_yaxes(title_text='Substitution Frequency', range=[-5.5,3.5],
                     tickmode = 'array',
                     tickvals = [x for x in range(-6,4,1)],
                     ticktext = [str(10**x) for x in range(-6,4,1)]
                     )
    return fig



#%%
@app.callback(
    Output('volcano','figure'),
    Input('volcanoType','value'),
    Input('volcanoNumerator','value'),
    Input('volcanoDenominator','value'),
    )
def updateVolcanoPlot(volcanoType,volcanoNumerator,volcanoDenominator):
    #Parse relevant data
    dff = pd.concat([bySample.get_group(name) for name in [volcanoNumerator,volcanoDenominator]])
    bySampleMetric = dff.groupby(by=['Sample',volcanoType])['Normalized'].apply(list)
    bySampleMetric = bySampleMetric.reset_index().pivot_table(columns='Sample',index=volcanoType,values='Normalized',
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
    maskSignificant = bySampleMetric['P'] <= benpv
    maskSignificant = (maskSignificant & (np.abs(bySampleMetric['FC']) > 1))
    bySampleMetric['logP'] = np.log10(bySampleMetric['P'])*-1
    sigmin = bySampleMetric['logP'].min()
    sigmax = bySampleMetric['logP'].max()
    fcmin = bySampleMetric['FC'].min()*missingoffset
    fcmax = bySampleMetric['FC'].max()*missingoffset

    fig = go.Figure()
    #Add significantly enriched values
    fig.add_trace(go.Scatter(x=bySampleMetric.loc[maskSignificant & maskBoth,'FC'],
                             y=bySampleMetric.loc[maskSignificant & maskBoth,'logP'],
                             mode='markers',marker_color='rgba(255,0,0)',
                            name='Significant',
                            hoverinfo='text',
                            text=bySampleMetric.loc[maskSignificant].index)
                  )
    #Add values only in one sample
    #random x is to add jitter
    fig.add_trace(go.Scatter(x=[fcmin*(1+random.randrange(1,10)/100) for x in range(len(bySampleMetric.loc[maskNaNumerator]))],
                             y=np.arange(sigmin,sigmax,sigmax/len(bySampleMetric.loc[maskNaNumerator])),
                             mode='markers',marker_color='rgba(255,0,0)',
                            name='Unique to ' + volcanoDenominator,
                            hoverinfo='text',
                            text=bySampleMetric.loc[maskNaNumerator].index)
                  )
    fig.add_trace(go.Scatter(x=[fcmax*(1+random.randrange(1,10)/100) for x in range(len(bySampleMetric.loc[maskNaDenominator]))],
                             y=np.arange(sigmin,sigmax,sigmax/len(bySampleMetric.loc[maskNaDenominator])),
                             mode='markers',marker_color='rgba(255,0,0)',
                            name='Unique to '+volcanoNumerator,
                            hoverinfo='text',
                            text=bySampleMetric.loc[maskNaDenominator].index)
                  )
    #Add insignificant values
    fig.add_trace(go.Scatter(x=bySampleMetric.loc[~maskSignificant & maskBoth,'FC'],
                             y=bySampleMetric.loc[~maskSignificant & maskBoth,'logP'],
                             mode='markers',marker_color='rgba(255,0,0)',
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
    fig.add_trace(go.Scatter(x=[-1,-1],y=[0,sigmax*1.2],
                             mode='lines',marker_color='rgba(130,130,130)',
                            showlegend=False,
                            line=dict(width=1,dash='dot'),
                            name='FC Neg line')
                  )
    #FC upper threshold
    fig.add_trace(go.Scatter(x=[1,1],y=[0,sigmax*1.2],
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
    fig.add_annotation(text="Higher "+'Normalized'+" in "+ volcanoNumerator,
                xref="paper", yref="paper",
                x=.85, y=-0.12,showarrow=False)
    fig.add_annotation(text="Higher "+'Normalized'+" in "+ volcanoDenominator,
                xref="paper", yref="paper",
                x=.15, y=-0.12,showarrow=False)
    
    fig.update_xaxes(title_text='log2(fold change)', range=[fcmin*1.2,fcmax*1.2],
                     
                     )
    fig.update_yaxes(title_text='P', range=[0,round(sigmax)+1],
                     tickmode = 'array',
                     tickvals = [x for x in range(0,round(sigmax)+1)],
                     ticktext = [str(10**-x) for x in range(0,round(sigmax)+1)]
                     )
    #dit
    return fig

#%%
if __name__ == '__main__':
    app.run_server(debug=True)