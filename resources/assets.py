#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:49:17 2021
Common resources
@author: taylorlundgren
"""
import os
import pandas as pd
import datetime
import numpy as np
package_dir = os.path.dirname(os.path.split(os.path.abspath(__file__))[0])
file_path = os.path.join(package_dir, 'temp', 'data_directory.txt')

#Variables you must change for certain scripts
Rscript = r"C:\Program Files\R\R-4.3.1\bin\Rscript.exe"     #For differential analysis with limma
#"C:\Program Files\R\R-4.3.2\bin\Rscript.exe" Fanta
#C:\Program Files\R\R-4.3.2\bin\Rscript.exe" laptop
#"C:\Program Files\R\R-4.3.1\bin\Rscript.exe" Celcius

#Paired DNA and protein fasta files 
    #MG1655
#path_dna_fasta = os.path.join(package_dir,r"resources\MG1655_dna_U00096_3_2024_04.fasta")
#path_protein_fasta = os.path.join(package_dir,r"resources\MG1655_protein_U00096_3_2024_04.fasta")
    #Xac
#path_dna_fasta = os.path.join(package_dir,r"resources\MG1655_dna_Xac_common_2024_05.fasta")
#path_protein_fasta = os.path.join(package_dir,r"resources\MG1655_protein_Xac_common_2024_05.fasta")
    #Yeast 
path_dna_fasta = os.path.join(package_dir,r'resources\W303_JRIU00000000_SGD_cds_trimmed2.fsa')
path_protein_fasta = os.path.join(package_dir,r'resources\W303_JRIU00000000_SGD_pep.fsa')


#For two-organism experiment
dbused = 'ECOLI'
targetlist = pd.read_csv(os.path.join(package_dir,'resources','SSP_ECOLI_to_SALTY.csv'))
targetlistrev = pd.read_csv(os.path.join(package_dir,'resources','SSP_SALTY_to_ECOLI.csv'))


#Directory management
def changeDir(newdir):
    with open(file_path, "w") as f:
        f.write(newdir)
    inputdir = newdir
    outputdir = os.path.join(newdir,'substitutuionAnnotation_output')


with open(file_path, "r+") as f:
    inputdir = f.read().strip()
    if not inputdir:
        inputdir = input('Please enter the directory with the data:')
        f.write(inputdir)
        

outputdir= os.path.join(inputdir,'substitutionAnnotation_output')
if not os.path.exists(outputdir):
    try:
        os.mkdir(outputdir)
    except:
        inputdir = input(inputdir +' is not a valid directory. Please try again.')
        
        try:
            os.mkdir(os.path.join(inputdir,'substitutionAnnotation_output'))
            changeDir(inputdir)
            outputdir= os.path.join(inputdir,'substitutionAnnotation_output')

        except FileExistsError:
            pass
        except:
            raise Exception('Cannot make that directory.')


#Commonly used lists and dictionaries

dfdm = pd.read_csv(os.path.join(package_dir,'resources','dfdm.csv'),index_col=0)
dfdm_varcys = pd.read_csv(os.path.join(package_dir,'resources','dfdm_varCysmod.csv'),index_col=0)
dfdangermods = pd.read_csv(os.path.join(package_dir,'resources','dangermods2.csv'))


dict_321 = {"Gly": 'G', 
            "Ala" : 'A', 
            "Ser" : 'S', 
            "Pro" : 'P', 
            "Val" : 'V', 
            "Thr" : 'T', 
            "Ile" : 'I',
            "Leu" : "L",
            "Asn" : 'N', 
            "Asp" : 'D', 
            "Gln" : 'Q', 
            "Lys" : 'K', 
            "Glu" : 'E', 
            "Met" : 'M', 
            "His" : 'H',
            "Phe" : 'F', 
            "Arg" : 'R', 
            "Cys" : 'C', 
            "Tyr" : 'Y',
            "Trp" : 'W',
            "Cis" : 'C',
            "Stp" : 'B'
            }

dict_123 = {'G':"Gly",
            'A':"Ala",
            'S':'Ser',
            'P':'Pro',
            'V':'Val',
            'T':'Thr',
            'I':'Ile',
            'L':'Leu',
            'N':'Asn',
            'D':'Asp',
            'Q':'Gln',
            'K':'Lys',
            'E':'Glu',
            'M':'Met',
            'H':'His',
            'F':'Phe',
            'R':'Arg',
            'C':'Cys',
            'Y':'Tyr',
            'W':'Trp',
            'B':'Stp'}

#Use a dictionary which assigns I and L to 'Xle' to reflect isobaric masses
dict_123ambiguous = {'G': 'Gly',
 'A': 'Ala',
 'S': 'Ser',
 'P': 'Pro',
 'V': 'Val',
 'T': 'Thr',
 'L': 'Xle',
 'I': 'Xle',
 'N': 'Asn',
 'D': 'Asp',
 'Q': 'Gln',
 'K': 'Lys',
 'E': 'Glu',
 'M': 'Met',
 'H': 'His',
 'F': 'Phe',
 'R': 'Arg',
 'C': 'Cys',
 'Y': 'Tyr',
 'W': 'Trp',
 'B':'Stp'}

dict_321ambiguous = {"Gly": 'G', 
            "Ala" : 'A', 
            "Ser" : 'S', 
            "Pro" : 'P', 
            "Val" : 'V', 
            "Thr" : 'T', 
            "Xle" : "L",
            "Asn" : 'N', 
            "Asp" : 'D', 
            "Gln" : 'Q', 
            "Lys" : 'K', 
            "Glu" : 'E', 
            "Met" : 'M', 
            "His" : 'H',
            "Phe" : 'F', 
            "Arg" : 'R', 
            "Cys" : 'C', 
            "Tyr" : 'Y',
            "Trp" : 'W',
            "CysCam" : 'C',
            "CysCamOx" : 'C',
            "Stp" : '*'
            }

def timestamp(text):
    print(text)
    now = datetime.datetime.now()
    print(str(now))
     

aminoacids_ambiguous = ['Ala', 'Arg','Asn','Asp','Cys','Glu','Gln','Gly',
              'His','Xle','Lys','Met','Phe','Pro','Ser',
              'Thr','Trp','Tyr','Val']

aminoacids = ['Ala', 'Arg','Asn','Asp','Cys','Glu','Gln','Gly',
              'His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
              'Thr','Trp','Tyr','Val']

#List of single letter AAs in alphabetical order, also in order of dna_codons below
aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

dna_codons = [
    ['GCT', 'GCC', 'GCA', 'GCG'],  # A
    ['TGT', 'TGC'],                # C
    ['GAT', 'GAC'],                # D
    ['GAA', 'GAG'],                # E
    ['TTT', 'TTC'],                # F
    ['GGT', 'GGC', 'GGA', 'GGG'],  # G
    ['CAT', 'CAC'],                # H
    ['ATT', 'ATC', 'ATA'],        # I
    ['AAA', 'AAG'],                # K
    ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # L
    ['ATG'],                       # M
    ['AAT', 'AAC'],                # N
    ['CCT', 'CCC', 'CCA', 'CCG'],  # P
    ['CAA', 'CAG'],                # Q
    ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # R
    ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # S
    ['ACT', 'ACC', 'ACA', 'ACG'],  # T
    ['GTT', 'GTC', 'GTA', 'GTG'],  # V
    ['TGG'],                       # W
    ['TAT', 'TAC']                 # Y
]

codon_table_inverse = pd.Series([
    ['GCT', 'GCC', 'GCA', 'GCG'],
    ['TGT', 'TGC'],                
    ['GAT', 'GAC'],               
    ['GAA', 'GAG'],               
    ['TTT', 'TTC'],                
    ['GGT', 'GGC', 'GGA', 'GGG'], 
    ['CAT', 'CAC'],               
    ['ATT', 'ATC', 'ATA'],      
    ['AAA', 'AAG'],                # K
    ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # L
    ['ATG'],                       # M
    ['AAT', 'AAC'],                # N
    ['CCT', 'CCC', 'CCA', 'CCG'],  # P
    ['CAA', 'CAG'],                # Q
    ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # R
    ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # S
    ['ACT', 'ACC', 'ACA', 'ACG'],  # T
    ['GTT', 'GTC', 'GTA', 'GTG'],  # V
    ['TGG'],                       # W
    ['TAT', 'TAC']                 # Y
], index = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])

codons = ['GCT', 'GCC', 'GCA', 'GCG','TGT', 'TGC','GAT', 'GAC', 'GAA',
 'GAG', 'TTT', 'TTC','GGT', 'GGC', 'GGA', 'GGG', 'CAT', 'CAC', 'ATT',
  'ATC', 'ATA', 'AAA', 'AAG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
  'ATG', 'AAT', 'AAC','CCT', 'CCC', 'CCA', 'CCG', 'TAA','TAG','TGA', 
    'CAA', 'CAG','CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',   
    'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC','ACT', 'ACC',
'ACA', 'ACG','GTT', 'GTC', 'GTA', 'GTG', 'TGG', 'TAT', 'TAC']

def clearDir():
    """Clears the .txt file used to store the working data directoty"""
    with open(file_path, "w") as f:
        f.write('')
        
def log10CV(series):
    """
    A function to return the %CV for log10 transformed data
    """
    return np.sqrt(10**(np.log(10)*np.std(series)**2)-1)*100

#Legacy code support fasta references

#Default is the MG1655 protein fasta
path_to_fasta = os.path.join(package_dir,'resources',r'2022-03-26-decoys-contam-UP_2022_03_25_EcoliK12.fasta.fas')
#path_to_fasta = os.path.join(package_dir,'resources',r'2021-11-mSmegMC2155.fasta')

#Default is to the MG1655 annotated genome GenBank
path_to_GenBank = os.path.join(package_dir,'resources',r'2023-10-20-GenBank-MG1655-Coding-DNA.gb')
#path_to_GenBank = os.path.join(package_dir,'resources',r'2023-10-26-smegmc2155.gb')
