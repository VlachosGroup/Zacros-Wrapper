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
RunPath = 'C:/Users/mpnun/Desktop/test_rep/1/'
ProductSpecies = 'CO2'
exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
''' -------------------------------------------- '''

''' Set up data '''
y = KMC_Run()
y.data.Path = RunPath
y.data.ReadAllInput()
y.exe_file = exe_file

#y.Run_sim()
y.data.ReadAllOutput()

''' Analyze '''
#TOF = y.ComputeTOF(ProductSpecies)
y.PlotSurfSpecVsTime(save = True)
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()
#y.LatticeMovie()