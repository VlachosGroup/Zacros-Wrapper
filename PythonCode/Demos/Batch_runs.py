# -*- coding: utf-8 -*-
"""
Created on Sun Oct 02 13:22:25 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling
from AnalyzeData import AnalyzeData
from KMCrun import KMCrun

os.system('cls')

# Set all directories
exe_path = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/'
KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB_nonstiff/'
RunPath = 'C:/Users/mpnun/Desktop/rescale_test/'

# Read in data
KMC_template = KMCrun()
KMC_template.data.Path = KMC_source
KMC_template.data.ReadAllInput()
KMC_template.data.Path = RunPath

KMC_batch = 