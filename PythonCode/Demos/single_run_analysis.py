# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from KMCrun import KMCrun


################################################################

os.system('cls')

#RunPath = 'C:/Users/mpnun/Desktop/test/'
#RunPath = 'C:/Users/mpnun/Desktop/WGS_run/'
#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/AtoB/0111/'
RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/AtoB/0001/'

''' Single run '''
y = KMCrun()
y.data.Path = RunPath
y.data.ReadAllInput()
print y.data.Conditions

#y.data.ReadAllOutput()
#TOF = y.ComputeTOF('B')
#y.PlotSurfSpecVsTime()
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()
#y.PlotWVsTime()
#y.PlotPropsVsTime()
#y.PlotIntPropsVsTime()