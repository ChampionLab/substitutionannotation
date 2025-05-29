# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:08:53 2022
A script to create an interactive plotly charts reflecting the abundance of substitution distribution
with interactive filtering
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



#%%
#Get data to plot
dfall = pd.read_csv(os.path.join(r'C:\Users\taylo\TEmp','AllPSMsAndFilters.csv'))
dfall = dfall[~dfall['Protein'].str.contains('cont|rev_')]                  #Remove contaminants and decoys
dfall['Intensity'].replace(0,np.nan,inplace=True)
psmSubs = dfall[dfall['Is Sub']]
psmSubs = psmSubs[~psmSubs['Is Danger']]                                    #Remove other PTMs

psmBase = dfall[dfall['Is Base']]

#Lists to choose from
annotations = ['Protein','Origin','Destination']
second_dimension = ['Ion Mobility shift','Retention shift','PeptideProphet Probability','Base Intensity','CV','Base CV','Ratio CV']
samples = psmSubs['Sample'].unique()

#Check/make directory for outputs
output = os.path.join(sA.outputdir,'dynamicVisualizations')
if not os.path.exists(output):
    os.mkdir(output)

#%%
#Start Dash

plotly.io.renderers.default='browser'
app = Dash(__name__,suppress_callback_exceptions=True)

#%%
#Layout of page

app.layout = html.Div([
    html.Div(children=[
        html.Div(children=[
            html.H1("Filters"),

            html.Label('Samples'),
            dcc.Checklist(id='sample_list',options=samples, value=samples),
            
            html.Label("Minimum intensity"),
            dcc.Input(id='min_intensity', type='number', value=0),

            html.Label("Minimum number of replicates"),
            dcc.Input(id='min_replicates', type='number', value=0),

            html.Label("Minimum number of samples"),
            dcc.Input(id='min_samples', type='number', value=0),
            
            html.Label("Minimum base intensity"),
            dcc.Input(id='min_base_intensity', type='number', value=0),
            
            html.Label("Maximum CV"),
            dcc.Input(id='max_cv', type='number', value=2),
            ]),
        html.Div(children=[
            html.Label("Maximum Base CV"),
            dcc.Input(id='max_base_cv', type='number', value=2),
            
            html.Label("Maximum ratio CV"),
            dcc.Input(id='max_ratio_cv', type='number', value=2),
            
            html.Label("Minimum ΔRT"),
            dcc.Input(id='min_delta_rt', type='number', value=0),
            
            html.Label("Minimum ΔIM"),
            dcc.Input(id='min_delta_im', type='number', value=0),
            
            html.Label("Minimum Probability"),
            dcc.Input(id='min_probability', type='number', value=0),
            
            html.Label("Protein"),
            dcc.Input(id='protein', type='text', value=''),
            
            html.Label("Origin"),
            dcc.Input(id='origin', type='text', value=''),
            
            html.Label("Destination"),
            dcc.Input(id='destination', type='text', value=''),

            html.Label('Hoverinfo'),
            dcc.Checklist(options=annotations, id='annotation_setting'),

            html.Label('Second dimension of data'),
            dcc.Dropdown(options=second_dimension, value=second_dimension[0], id='second_dimension'),
            ]),
        ], style={'display': 'flex', 'flexDirection': 'row'}),
        
    html.Div(children=[
        html.Hr(),  # Add a horizontal line
        
        html.H1("Plots"),
        
        html.Div([
            html.Div([
                html.H3("Violin Plot - Intensity"),
                dcc.Graph(id='violin_intensity')
            ], className="six columns"),
            
            html.Div([
                html.H3("Violin Plot - Ratio"),
                dcc.Graph(id='violin_ratio')
            ], className="six columns"),
            
            html.Div([
                html.H3("Density Plot"),
                dcc.Graph(id='density_plot')
            ], className="twelve columns")
        ], className="row")
    ]),
])


#%%

@app.callback(
    Output('violin_intensity', 'figure'),
    Output('violin_ratio', 'figure'),
    Output('density_plot', 'figure'),

    Input('min_intensity', 'value'),
    Input('min_base_intensity', 'value'),
    Input('max_cv', 'value'),
    Input('max_base_cv', 'value'),
    Input('max_ratio_cv', 'value'),
    Input('min_delta_rt', 'value'),
    Input('min_delta_im', 'value'),
    Input('min_probability', 'value'),
    Input('protein', 'value'),
    Input('origin', 'value'),
    Input('destination', 'value'),
    Input('annotation_setting', 'value'),
    Input('second_dimension', 'value'),
    Input('sample_list','value'),
    Input('min_replicates','value'),
    Input('min_samples','value') 
)
def filter_and_plot(min_intensity, min_base_intensity, max_cv, max_base_cv, max_ratio_cv, min_delta_rt,
                             min_delta_im, min_probability, protein, origin, destination, annotation_setting, second_dimension,
                             sample_list,min_replicates,min_samples):
    
    #Filter data
    psmSubsf = psmSubs[psmSubs['Sample'].isin(sample_list)]
    psmBasef = psmBase[psmBase['Sample'].isin(sample_list)]
    psmSubsf = psmSubsf[psmSubsf['Protein'].str.contains(protein)]
    psmBasef = psmBasef[psmBasef['Protein'].str.contains(protein)]
    psmSubsf = psmSubsf[psmSubsf['Origin'].str.contains(origin)]
    psmSubsf = psmSubsf[psmSubsf['Destination'].str.contains(destination)]
    psmSubsf = psmSubsf[psmSubsf['Intensity'] >= min_intensity]
    psmBasef = psmBasef[psmBasef['Intensity'] >= min_base_intensity]
    psmSubsf = psmSubsf[psmSubsf['PeptideProphet Probability'] >= min_probability]
    psmBasef = psmBasef[psmBasef['PeptideProphet Probability'] >= min_probability]

    #Aggregate replicate ata
    byRepSubs = psmSubsf.groupby(by=['Modified Peptide','Sample','Replicate'])[['Intensity','Retention','Ion Mobility']].mean()
    byRepSubs[['Peptide','Protein','Origin','Destination']] = psmSubsf.groupby(by=['Modified Peptide','Sample','Replicate'])[['Peptide','Protein','Origin','Destination']].apply(lambda x: x.iloc[0,:])   #Take first value
    byRepSubs.reset_index(inplace=True)
    byRepBase = psmBasef.groupby(by=['Peptide','Sample','Replicate'])[['Intensity','Retention','Ion Mobility']].mean().reset_index()

    byRepBase.rename(columns={'Intensity':'Base Intensity','Retention':'Base Retention','Ion Mobility':'Base Ion Mobility'},inplace=True)

    byRep = byRepSubs.merge(byRepBase,how='left',on=['Peptide','Sample','Replicate'])
    byRep['Retention shift'] = byRep['Retention'] - byRep['Base Retention']
    byRep['Ion Mobility shift'] = byRep['Ion Mobility'] - byRep['Base Ion Mobility']

    #Apply replicate filters
    byRep = byRep[np.abs(byRep['Retention shift']) >= min_delta_rt]
    byRep = byRep[np.abs(byRep['Ion Mobility shift']) >= min_delta_im]

    #Aggregate sample data
    a = byRep.groupby(by=['Modified Peptide','Sample'])[['Intensity','Retention','Ion Mobility','Base Intensity','Retention shift','Ion Mobility shift']].mean()
    b = byRep.groupby(by=['Modified Peptide','Sample'])[['Peptide','Protein','Origin','Destination']].max()
    c = byRep.groupby(by=['Modified Peptide','Sample'])[['Intensity','Base Intensity']].apply(lambda x: x.std()/x.mean())
    c.columns = ['CV','Base CV']
    d = byRep.groupby(by=['Modified Peptide','Sample'])['Intensity'].count().rename('N Replicates')
    bySample = pd.concat([a,b,c,d], axis=1)
    bySample = bySample.reset_index()
    bySample['Ratio'] = bySample['Intensity']/bySample['Base Intensity']
    bySample['logRatio'] = np.log10(bySample['Ratio'])
    bySample['Ratio CV'] = np.sqrt(bySample['CV']**2 + bySample['Base CV']**2)

    #Apply sample filters
    bySample = bySample[bySample['CV'] < max_cv]
    bySample = bySample[bySample['Base CV'] < max_base_cv]
    bySample = bySample[bySample['Ratio CV'] < max_ratio_cv]
    bySample = bySample[bySample['N Replicates'] >= min_replicates]

    #Figure out n samples
    nsamples = bySample.groupby(by=['Modified Peptide'])['Intensity'].count().rename('N Samples').reset_index()
    
    bySample = bySample.merge(nsamples,how='left',on='Modified Peptide')
    bySample = bySample[bySample['N Samples'] > min_samples]
    allSamples = bySample.copy()
    bySample = bySample.groupby(by=['Sample'])
    #Plot data

    #Intesnity scatter
    intensity_scatter = px.scatter(allSamples,x='Base Intensity',y='Intensity',
                        facet_col='Sample',log_x=True,log_y=True,hover_data=annotation_setting)

    #Ratio vioin

    violin_ratio = px.violin(allSamples,x='Sample',y='logRatio',hover_data=annotation_setting,
        box=True,points='all')
    violin_ratio.update_layout(violingap=0,showlegend=False)

    #Density plot
    # Create density heatmap
    density_plot = px.scatter(allSamples, x=second_dimension, y="logRatio", facet_col='Sample',hover_data=annotation_setting)

    return intensity_scatter, violin_ratio, density_plot




#%%
if __name__ == '__main__':
    app.run(debug=True)