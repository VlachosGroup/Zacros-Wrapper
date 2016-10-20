# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from KMC_Run import KMC_Run

os.system('cls')

''' ------------ User input section ------------ '''
RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/lattice_test/'
ProductSpecies = 'B'
''' -------------------------------------------- '''

''' Set up data '''
y = KMC_Run()
y.data.Path = RunPath
y.data.ReadAllInput()
y.data.ReadAllOutput()

''' Analyze '''
TOF = y.ComputeTOF(ProductSpecies)
y.PlotSurfSpecVsTime()
y.PlotGasSpecVsTime()
y.PlotElemStepFreqs()
y.LatticeMovie()