# -*- coding: utf-8 -*-
"""
Created on Fri Nov 04 19:02:14 2016

@author: mpnun
"""

import os
import sys
import copy

sys.path.insert(0, '../KMCsim')
from KMC_Run import KMC_Run

os.system('cls')

''' ------------ User input section ------------ '''
RunPath = 'C:/Users/mpnun/Desktop/test_cont/1/'
ProductSpecies = 'B'
exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
''' -------------------------------------------- '''

''' Set up data '''
y = KMC_Run()
y.Path = RunPath
y.exe_file = exe_file
y.ReadAllInput()

#y.Run_sim()
y.ReadAllOutput()
#y.PlotSurfSpecVsTime(save = False)

z = copy.deepcopy(y)
z.Path = 'C:/Users/mpnun/Desktop/test_cont/2/'
z.StateInput['Type'] = 'history'
z.StateInput['Struct'] = y.History[-1]
#z.WriteAllInput()
#z.Run_sim()
z.ReadAllOutput()
#z.PlotSurfSpecVsTime(save = False)

sw = KMC_Run.time_sandwich([y,z])


#sw.PlotSurfSpecVsTime(save = False)

y.PlotPropsVsTime()


''' Analyze '''
#TOF = y.ComputeTOF(ProductSpecies)
#y.PlotSurfSpecVsTime(save = False)
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()
#y.LatticeMovie()