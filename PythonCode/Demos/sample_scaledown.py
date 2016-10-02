# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

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
KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB/'
RunPath = 'C:/Users/mpnun/Desktop/rescale_test/'

''' Single run '''
y = KMCrun()
y.data.Path = KMC_source
y.data.ReadAllInput()
y.data.Path = RunPath

#y.data.ReadAllOutput()
z = RateRescaling()
z.KMC_system = y

#delta_sdf = z.ProcessStepFreqs()
#print delta_sdf

z.PerformScaledown(exe_path)
z.KMC_system.PlotElemStepFreqs()
z.PlotStiffnessReduction()
z.WriteRescaling_output()
z.KMC_system.CheckSteadyState('B', show_graph = True)
print 'Final KMC time: ' + str(z.KMC_system.data.Specnum['t'][-1])