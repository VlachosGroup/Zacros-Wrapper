# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

from KMCrun import KMCrun
from OutputData import OutputData
from AnalyzeData import AnalyzeData
import os
import numpy as np
import matplotlib.pyplot as plt

#import sys
#import subprocess
#import pickle
#import matplotlib.pyplot as plt

################################################################

os.system('cls')

#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/AtoB/0111/'
#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/WGS/111/'
#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/LocalRun/Run/'
RunPath1 = 'C:/Users/mpnun/Desktop/fldr1/'
RunPath2 = 'C:/Users/mpnun/Desktop/fldr2/'

''' Single run '''
y = KMCrun()
y.data.Path = RunPath1
y.data.ReadAllOutput()

y.data.Path = RunPath2
y.data.WriteAllInput()

y.PlotWVsTime()

#print 'Mode'
#print y.output.input.StiffnessRecondition['Mode']
#
#print 'Scaledown factor'
#print y.output.input.StiffnessRecondition['SDF']
#
##print y.output.input.Reactions['Input']
#print y.output.input.Reactions['nrxns']

#TOF = y.ComputeTOF('B')
#y.PlotSurfSpecVsTime()
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()
#y.PlotWVsTime()
#y.PlotPropsVsTime()
#y.PlotIntPropsVsTime()