# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from KMC_batch import KMC_batch

################################################################

os.system('cls')

BatchPath = 'C:/Users/mpnun/Desktop/analyzethese/'

# Batch of runs ----------------
x = KMC_batch()
x.ParentFolder = BatchPath
x.ReadMultipleRuns()

# Trajectory average
x.AverageRuns()
x.runAvg.PlotSurfSpecVsTime()
x.runAvg.PlotGasSpecVsTime()
x.runAvg.PlotElemStepFreqs()
#x.runAvg.CheckSteadyState('CO2', show_graph = True)

# Rate and sensitivities
x.ComputeStats('CO2')
#x.WvarCheck() 
x.PlotSensitivities()
x.WriteSA_output(BatchPath)