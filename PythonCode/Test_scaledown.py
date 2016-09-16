# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

from KMCrun import KMCrun
from AnalyzeData import AnalyzeData
import os
import numpy as np
import matplotlib.pyplot as plt
from RateRescaling import RateRescaling

#import sys
#import subprocess
#import pickle
#import matplotlib.pyplot as plt

################################################################

os.system('cls')

#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/AtoB/0111/'
#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/JobBuilds/WGS/111/'
#RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/LocalRun/Run/'
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

z.PerformScaledown()
z.PlotStiffnessReduction()