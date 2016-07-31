# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

#from KMCrun import KMCrun
from OutputData import OutputData

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
y = OutputData()
y.Path = RunPath
y.ReadAllOutput()

print y.Binary['W_sen_anal']