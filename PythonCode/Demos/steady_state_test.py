# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:32:59 2016

@author: mpnun
"""

from KMCrun import KMCrun
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

BatchPath = 'C:/Users/mpnun/Desktop/not_ss/'

''' Group of runs '''

# Test batch of runs ----------------
x = AnalyzeData()
x.ReadMultipleRuns(BatchPath)
x.AverageRuns()


#x.runAvg.PlotSurfSpecVsTime()
#x.runAvg.PlotGasSpecVsTime()
#x.runAvg.PlotElemStepFreqs()

x.ComputeStats('B')
is_ss = x.runAvg.CheckSteadyState('B')
print is_ss
#x.PlotSensitivities()
#x.WriteSA_output(BatchPath)