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

package_dir = os.path.dirname(os.path.split(os.path.abspath(__file__))[0])
file_path = os.path.join(package_dir, 'temp', 'data_directory.txt')


with open(file_path, "r+") as f:
    inputdir = f.read().strip()
    if not inputdir:
        inputdir = input('Please enter the directory with the data:')
        f.write(inputdir)
        

outputdir= os.path.join(inputdir,'AASoutput')
if not os.path.exists(outputdir):
    os.mkdir(outputdir)

#Commonly used lists

dfdm = pd.read_csv(os.path.join(package_dir,'resources','dfdm.csv'),index_col=0)
dfdangermods = pd.read_csv(os.path.join(package_dir,'resources','dangermods.csv'))


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
            "Cis" : 'C'
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
            'W':'Trp'}

#Use a dictionary which assigns I and L to 'Leu/Ile' to reflect isobaric masses
dict_123ambiguous = {'G': 'Gly',
 'A': 'Ala',
 'S': 'Ser',
 'P': 'Pro',
 'V': 'Val',
 'T': 'Thr',
 'L': 'Leu/Ile',
 'I': 'Leu/Ile',
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
 'W': 'Trp'}

def timestamp(text):
    print(text)
    now = datetime.datetime.now()
    print(str(now))
     

aminoacids = ['Ala', 'Arg','Asn','Asp','Cys','Glu','Gln','Gly',
              'His',u'Leu/Ile','Lys','Met','Phe','Pro','Ser',
              'Thr','Trp','Tyr','Val']


