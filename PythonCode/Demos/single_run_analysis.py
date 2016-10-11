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

RunPath = 'C:/Users/mpnun/Desktop/analyzethese/1/'

''' Single run '''
y = KMCrun()
y.data.Path = RunPath
y.data.ReadAllInput()

#print len(y.data.Reactions['Input'][0]['variant'])

#for x in y.data.Reactions['Input'][0]:
#    print (x) + ': ' + str(y.data.Reactions['Input'][0][x])
#    print y.data.Reactions['Input'][0][x]
#    for y in y.data.Reactions['Input'][0][x]:
#        print (y,':', y.data.Reactions['Input'][0][x][y])

print y.data.Reactions['names']

#print y.data.Reactions['Input'][0]['Name']

#y.data.ReadAllOutput()
#TOF = y.ComputeTOF('B')
#y.PlotSurfSpecVsTime()
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()
#y.PlotWVsTime()
#y.PlotPropsVsTime()
#y.PlotIntPropsVsTime()