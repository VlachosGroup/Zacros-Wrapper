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

RunPath = 'C:/Users/mpnun/Desktop/lattice_test/'

''' Single run '''
y = KMCrun()
y.data.Path = RunPath
y.data.ReadAllInput()
y.data.ReadAllOutput()

y.LatticeMovie()

#print y.data.History

#TOF = y.ComputeTOF('B')
#y.PlotSurfSpecVsTime()
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()
#y.PlotWVsTime()
#y.PlotPropsVsTime()
#y.PlotIntPropsVsTime()