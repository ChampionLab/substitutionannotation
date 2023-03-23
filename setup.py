# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 16:49:31 2023

@author: tjl
"""

from setuptools import setup, find_packages

setup(
    name='proteomicAnalysis',
    version='0.001',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'seaborn'
    ],
)
