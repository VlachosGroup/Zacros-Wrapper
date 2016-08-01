# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

from KMCrun import KMCrun
from OutputData import OutputData
from AnalyzeData import AnalyzeData

#import sys
import os
#import numpy as np
#import subprocess
#import pickle
#import matplotlib.pyplot as plt

################################################################

os.system('cls')
RunPath = 'C:/Users/mpnun/Desktop/test/1/'

## Test output ----------------
y = KMCrun()
y.output.Path = RunPath
y.output.ReadAllOutput()
TOF = y.ComputeTOF('B')
print TOF

#
#print y.Binary['W_sen_anal']

## Test analysis ----------------
#x = AnalyzeData()
#x.ReadMultipleRuns(RunPath)
#x.runList[0].ComputeTOF('B')