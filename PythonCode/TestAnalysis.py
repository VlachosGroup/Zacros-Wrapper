# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

from BuildInputFiles import BuildInputFiles as BI
from KMCUtilities import KMCUtilities as KMCut
from ProcessOutput import ProcessOutput as PO
from ReadInputFiles import ReadInputFiles as RI
from ReadOutputFiles import ReadOutputFiles as RO
from ReconditionStiffness import ReduceStiffness as RS
from RunZacros import RunZacros
from GeneralUtilities import GeneralUtilities as ut
import sys
import os
import numpy as np
import subprocess
import pickle
import matplotlib.pyplot as plt

os.system('cls')
RunPath = 'C:/Users/mpnun/Desktop/test/'

#Cnd = RO().ReadRun(RunPath)
#ut().PrintDict(Cnd)

#CndList = RO().ReadMultipleRuns(RunPath)
#print CndList[0]['Conditions']

x = KMCut()
print type(x)
print x.Info