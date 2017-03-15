# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.append('/home/vlachos/mpnunez/Github/Zacros-Wrapper/PythonCode')
import Core as zw

''' ------------ User input section ------------ '''
RunPath = '/home/vlachos/wangyf/CO_oxidation/23'
#RunPath = '/home/vlachos/mpnunez/ZW_data/KMC_data/WGS/SteadyState/Iteration_7/63'
''' -------------------------------------------- '''

''' Set up data '''
y = zw.KMC_Run()
y.Path = RunPath
#y.ReadAllInput()
#y.exe_file = exe_file

#y.Run_sim()
y.ReadAllOutput()
y.Path = '/home/vlachos/mpnunez/test'

''' Analyze '''
#n_Pd = 48 * 48
n_Pd = 920
y.PlotSurfSpecVsTime(site_norm = n_Pd)
y.PlotGasSpecVsTime()
y.PlotElemStepFreqs(site_norm = n_Pd, time_norm = True)
y.PlotLattice()
y.LatticeMovie()