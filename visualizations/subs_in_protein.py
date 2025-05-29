
# -*- coding: utf-8 -*-
"""
Created on  ?? 2024
A script to visualize substitution data on the 3D structure in PyMOL.
Must follow FindSubs_xxx.py
@author: taylo
"""

#%%
from dash import Dash, html, dcc, Input, Output, ctx
import pymol
from pymol import cmd
import os
from substitutionannotation.resources import assets as sA
import requests
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import re
from plotly import subplots as psubplots
import plotly
import plotly.graph_objects as go
import scipy.stats as stats

#%% Import data, pregroup data

global uniprot_id, alphafold_pdb, allbyprotf, subsbyprotf, dff

try:
    df = pd.read_csv(os.path.join(sA.outputdir,'Filtered SSP Quant.csv'))
except FileNotFoundError:
    df = pd.read_csv(os.path.join(sA.outputdir,'SSP Quant.csv'))

samplelist = df['Sample'].unique().tolist() + ['All']
dfannotate = pd.read_csv(os.path.join(sA.outputdir,'SSP PSM.csv'))
proteinlist = dfannotate['Protein'].drop_duplicates()
dfannotate.drop(columns=['Protein','Peptide','Intensity'],inplace=True)
eftu = proteinlist[proteinlist.str.contains('EFTU')].values[0]
dfall = pd.read_csv(os.path.join(sA.outputdir,'AllPSMsAndFilters.csv'))
allbysample = dfall.groupby(by='Sample')
allbyprot = dfall.groupby(by='Protein')
basebysample = dfall[dfall['PTMs']=='NONE'].groupby(by='Sample')
basebyprot = dfall[dfall['PTMs']=='NONE'].groupby(by='Protein')

df = df.merge(dfannotate,how='left',on=['Modified Peptide','Sample'])
subsbysample = df.groupby(by='Sample')
subsbyprot = df.groupby(by='Protein')
proteins = df['Protein'].unique()


dictcolor = {'Frequency Weighted Mean':'Substitution Sum Fraction',
             'Frequency Harmonic Mean':'Substitution Ratio Harmonic Mean',
             'Peptide Coverage':'Peptide Coverage'}

#%%
# PyMOL visualization functions

#Import protein
def pymol_protein_visualization(alphafold_pdb):
    """
    Download the PDB 3d file for uniprot, start pymol, and color the protein.
    alphafold_pdb = the PDB identifier to access api and get pdb file
    color_vector = Optional values used to color different residues.
    """
    name = re.findall('F\-(.*)\-F',alphafold_pdb)[0]

    cmd.load(alphafold_pdb,object=name)               #Get the PDB
    cmd.show('cartoon', name)     
    cmd.color('white', name)      #Default protein color
    cmd.zoom(name)
    
def pymol_color_residues(ser,alphafold_pdb):
    """
    Takes a series with aa position as index and numerical values.
    Converts numerical values based on a colormap, then colors the provided protein.
    """
    name = re.findall('F\-(.*)\-F',alphafold_pdb)[0]
    #Convert numerical values into colors
    cmap = plt.get_cmap('Wistia')
    cmap.set_bad('lightgrey')  #NaN color
    norm = mcolors.Normalize(vmin=ser.min(), vmax=ser.max())  

    colors = ser.apply(lambda x:cmap(norm(x)))   #Get (R,G,B,a)
    colors = colors.apply(rgba_to_rgb)

    #Color the protein
    for res in colors.index:
        cmd.color(colors[res],f'{name} and resi {res}')

def rgba_to_rgb(rgba):
    r,g,b,_ = rgba
    r = round(r*100)
    g = round(g*100)
    b = round(b*100)
    if r == 100:
        r = 99
    if g == 100:
        g = 99
    if b == 100:
        b = 99
    rgbstr = f'0x{r:02}{g:02}{b:02}'

    return rgbstr

def hex_to_int(hex_color_code): #Get integer representation of hex value
    r = int(hex_color_code[1:3], 16)
    g = int(hex_color_code[3:5], 16)
    b = int(hex_color_code[5:7], 16)
    int_representation = (r << 16) + (g << 8) + b
    return f'0x{int_representation}'


def get_protein(protein):

    uniprot_id = re.findall('\|(.*)\|',protein)[0]

    #Import data from Uniprot API
    response = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}")
    data = response.json()
    features = data.get("features", [])
    sequence = data.get("sequence")

    #Set up dataframe to hold protein information, with index as aa position (starting at 1)
    cols = ['Peptide Coverage','Substitution Coverage','Substitution Sum Fraction',
            'Sequence','Secondary Structure']
    dfprotein = pd.DataFrame(index=range(1,sequence['length']+1),columns=cols)

    #Get sequence
    dfprotein['Sequence'] = list(sequence['sequence'])

    #Get secondary structure
    for feature in features:
        if feature.get("category") == 'STRUCTURAL':
            dfprotein.loc[int(feature.get("begin")):int(feature.get("end")),'Secondary Structure'] = feature.get("type")

    #Get peptide coverage
    dfprotein['Peptide Coverage'] = False
    allpep = allbyprotf.get_group(protein)
    allpep = allpep.drop_duplicates(subset = ['Protein Start','Protein End'])
    for _,row in allpep.iterrows():
        dfprotein.loc[row['Protein Start']:row['Protein End'],'Peptide Coverage'] = True


    

    #Get substitution frequency by position
    dfprotein['Substitution Sum Fraction'] = get_sub_fraction(protein,protein_length = sequence['length']+1)

    #Get substitution ratio harmonic mean by position
    dfprotein['Substitution Ratio Harmonic Mean'] = subsbyprotf.get_group(protein).groupby(by=['Substitution Position'])['Ratio'].apply(stats.hmean)


    #Set alphafold_pdb 
    global alphafold_pdb 
    alphafold_pdb = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'

    return dfprotein

def get_sub_fraction(protein,protein_length):
    #Parse relevent protein data
    dfsubsf = subsbyprotf.get_group(protein)
    #Sum intensity of all substitutions that have the same base peptide and substitution position
    position_sumofsubs = dfsubsf.groupby(by=['Substitution Position','Peptide'])['Intensity'].sum()
    #Get the cognate intensity at each unique substitution position
    position_cognate_individual = dfsubsf.groupby(by=['Substitution Position','Peptide'])['Base Intensity'].max()
    #Perform weighted summing of substitution fraction at the same position
    weighted_sum = (position_sumofsubs/(position_sumofsubs+position_cognate_individual))*position_cognate_individual
    weighted_sum.name = 'Weighted Intensity'
    weighted_sum = weighted_sum.reset_index().groupby(by='Substitution Position')['Weighted Intensity'].sum()

    #Get total cognate signal at each position
    position_cognate_total = pd.Series(index=range(1,protein_length),data=[0]*(protein_length-1))
    basepeptides = basebyprotf.get_group(protein)
    for _,row in basepeptides.iterrows():
        position_cognate_total.loc[row['Protein Start']:row['Protein End']] += row['Intensity']

    sub_fraction = pd.Series(index=range(1,protein_length),data=[0]*(protein_length-1))
    sub_fraction = weighted_sum/position_cognate_total

    return sub_fraction

def mpl_to_plotly(cmap, pl_entries=11, rdigits=2):
    #Converts a matplotlib cmap to a plotly cmap
    # cmap - colormap 
    # pl_entries - int = number of Plotly colorscale entries
    # rdigits - int -=number of digits for rounding scale values
    scale = np.linspace(0, 1, pl_entries)
    colors = (cmap(scale)[:, :3]*255).astype(np.uint8)
    pl_colorscale = [[round(s, rdigits), f'rgb{tuple(color)}'] for s, color in zip(scale, colors)]
    return pl_colorscale

#%% Open/close pymol functions



#%% Dash interactive app

plotly.io.renderers.default='browser'
app = Dash(__name__,suppress_callback_exceptions=True)

app.layout = html.Div([
    html.H1("Protein Visualization"),
    dcc.Dropdown(options=samplelist,value='All',id='samples'),
    dcc.Dropdown(options=proteinlist,value=eftu,id='protein'),
    dcc.Graph(id='heatmap'),
    dcc.Dropdown(options=['Frequency Weighted Mean','Frequency Harmonic Mean','Peptide Coverage'],id='colorType'),
    html.Button('Launch PyMOL',id='openPymol'),
    html.Button('Color the 3D structure in PyMOL',id='colorbutton'),
    html.H3('',id='dummyout')
])

@app.callback(Output('colorbutton', 'value'),
              Input('samples','value'))
def sample_filter(samples):
    global allbyprotf, subsbyprotf, basebyprotf
    if samples == 'All':
        allbyprotf = allbyprot
        subsbyprotf = subsbyprot
        basebyprotf = basebyprot
    else:
        allbyprotf = allbysample.get_group(samples).groupby(by='Protein')
        subsbyprotf = subsbysample.get_group(samples).groupby(by='Protein')
        basebyprotf = basebysample.get_group(samples).groupby(by='Protein')



@app.callback(Output('heatmap','figure'),
              Input('protein','value'))
def update_heatmap(protein):
    df = get_protein(protein)
    df['Peptide Coverage'] = df['Peptide Coverage'].astype(int)
    dict_structure_int = {'HELIX':1,'STRAND':2,'TURN':3}
    df['Secondary Structure'] = df['Secondary Structure'].map(dict_structure_int).fillna(0)
    global dff
    dff = df.copy()
    #plot
    #Define colorbars
    wistia = plt.get_cmap('Wistia')
    colorgradient = mpl_to_plotly(wistia,250)
    colorstructure = [[0,'grey'],[.25,'grey'],
                      [.25,'green'],[.499,'green'],
                      [.499,'blue'],[.75,'blue'],
                      [.75,'orange'],[1,'orange']]
    colorcoverage = [[0,'grey'],[.5,'grey'],
                     [.51,'black'],[1,'black']]

    fig = psubplots.make_subplots(rows=4, cols=1,
                        row_heights=[0.2,0.2,0.2,0.2],
                        vertical_spacing = 0.005,
                        shared_xaxes=True,
                        )
    
    #Grey background for substitution frequency
    fig.add_trace(go.Heatmap(x=df.index,y=['Weighted Mean']*len(df),z=df['Substitution Sum Fraction'].fillna(0),colorscale=['grey','orange']),
            row=1,col=1)  
    #Substitution Frequency   
    fig.add_trace(go.Heatmap(x=df.index,y=['Weighted Mean']*len(df),z=np.log10(df['Substitution Sum Fraction']),colorscale=colorgradient,name='Ratio Weighted Mean'),
                row=1,col=1)

    #Grey background for harmonic frequency
    fig.add_trace(go.Heatmap(x=df.index,y=['Harmonic Mean']*len(df),z=df['Substitution Ratio Harmonic Mean'].fillna(0),colorscale=['grey','orange']),
            row=2,col=1)  
    #Harmonic Frequency 
    fig.add_trace(go.Heatmap(x=df.index,y=['Harmonic Mean']*len(df),z=df['Substitution Ratio Harmonic Mean'],colorscale=colorgradient,name='Ratio Harmonic Mean'),
                row=2,col=1)

    #Protein Coverage
    fig.add_trace(go.Heatmap(x=df.index,y=['Peptide Coverage']*len(df),z=df['Peptide Coverage'],colorscale=colorcoverage,name='Peptide Coverage'),
                row=3,col=1)

    #Protein Secondary Structure
    fig.add_trace(go.Heatmap(x=df.index,y=['Secondary Structure']*len(df),z=df['Secondary Structure'],colorscale=colorstructure,name='Secondary Structure'),
                row=4,col=1)

    #Annotations
    fig.update_traces(showscale=False)

    fig.update_traces(colorbar=dict(len=0.23,
                        orientation='v',
                        ticktext=['e-1','e-2','e-3','e-4','e-5'],
                        tickvals=[-1,-2,-3,-4,-5],
                        y=.875),
                    selector=dict(name='Ratio Weighted Mean'),
                    showscale=True)
    fig.update_traces(colorbar=dict(len=0.23,
                        orientation='v',
                        ticktext=['e-1','e-2','e-3','e-4','e-5'],
                        tickvals=[-1,-2,-3,-4,-5],
                        y=.625),
                    selector=dict(name='Ratio Harmonic Mean'),
                    showscale=True)
    fig.update_traces(colorbar=dict(len=0.23,
                        orientation='v',
                        tickmode='array',
                        ticktext=['Missing','Observed'],
                        tickvals=[.25,.75],
                        y=.375),
                    selector=dict(name='Peptide Coverage'),
                    showscale=True)
    fig.update_traces(colorbar=dict(len=0.23,
                        orientation='v',
                        tickmode='array',
                        ticktext=['Loop','Helix','Strand','Turn'],
                        tickvals=[0.37,1.12,1.87,2.62],
                        y=.125),
                    selector=dict(name='Secondary Structure'),
                    showscale=True)
    fig.update_xaxes(tickvals=df.index,ticktext=df['Sequence'])

    return fig

@app.callback(Output("openPymol","value"),
        Input("openPymol","n_clicks"))
def open_pymol(n_clicks):
    if "openPymol" == ctx.triggered_id:
        pymol.finish_launching()
    return 'PyMOL Launched'

@app.callback(Output("dummyout",'title'),
            Input("colorType","value"),
              Input('colorbutton','n_clicks'))
def open_pymol(colorType,n_clicks):
    if "colorbutton" == ctx.triggered_id:
        global dff, alphafold_pdb
        name = re.findall('F\-(.*)\-F',alphafold_pdb)[0]

        existing_objects = pymol.cmd.get_object_list()
        if not name in existing_objects:   #Don't open a fresh copy of the protein every time you color it
            pymol_protein_visualization(alphafold_pdb)

        pymol_color_residues(dff[dictcolor.get(colorType)],alphafold_pdb)

    return 'Ok'

if __name__ == "__main__":
    app.run_server(debug=True)
