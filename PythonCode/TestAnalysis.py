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
RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/WGS/111/'
BatchPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/WGS/'
#BatchPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/AtoB/'
#BatchPath = 'C:/Users/mpnun/Desktop/test/'

## Test single run ----------------
#y = KMCrun()
#y.output.Path = RunPath
#y.output.ReadAllOutput()
#TOF = y.ComputeTOF('CO2')
#y.PlotSurfSpecVsTime()
#y.PlotGasSpecVsTime()
#y.PlotElemStepFreqs()

# Test batch of runs ----------------
x = AnalyzeData()
x.ReadMultipleRuns(BatchPath)
x.AverageRuns()
#x.ComputeStats('B')
x.ComputeStats('CO2')
x.PlotSensitivities()