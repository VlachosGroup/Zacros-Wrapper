# -*- coding: utf-8 -*-
"""
Created on Tue Nov 08 15:52:33 2016

@author: mpnun
"""

import os
import sys

import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np

sys.path.insert(0, '../KMCsim')
from Replicates import Replicates

################## User input ##################################

BatchPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_scaledown2/Iteration_7/'
#BatchPath = 'C:/Users/mpnun/Desktop/test_rep/'
Product = 'B'
n_cores = 3

################################################################

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')

#    # Batch of runs ----------------
#    x = Replicates()
#    x.ParentFolder = BatchPath
#    x.n_procs = n_cores
#    x.ReadMultipleRuns(parallel = True)
#
#    # Trajectory average
#    x.AverageRuns()
#    
#    x.runAvg.CalcRateTraj(Product)
#    x.runAvg.PlotRateVsTime()
#    print x.runAvg.CheckSteadyState(Product)
#    print x.CheckAutocorrelation(Product)

    delt_list = []
    corr_list = []
    corr_err = []
    corr_expect_list = []
    for ind in range(2,5):
        BatchPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/ab/Iteration_' + str(ind) + '/'
        
        x = Replicates()
        x.ParentFolder = BatchPath
        x.n_procs = n_cores
        x.ReadMultipleRuns(parallel = True)
    
        # Trajectory average
        x.AverageRuns()
        
        x.runAvg.CalcRateTraj(Product)
#        print x.runAvg.CheckSteadyState(Product)
        corrs =  x.CheckAutocorrelation(Product)     
        
        delt = x.runAvg.Specnum['t'][-1]/2
        delt_list.append(delt)
        corr_list.append(corrs[0])
        corr_err.append(corrs[1])
        
 
    tau = 0.3186
    dtlin = np.linspace(0, 2, 1000)
    corr_expected = []
    for dt in dtlin:
        corr_expected.append(np.exp(-dt / tau)) 
    
    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 12    
    
    plt.figure()
#    plt.errorbar(range(2,8), corr_list, yerr = corr_err, markersize = 15, marker='o')
    plt.errorbar(delt_list, corr_list, yerr = corr_err, markersize = 15, marker='o')
    plt.plot(dtlin, corr_expected)
    plt.xticks(size=24)
    plt.yticks(size=24)
#    plt.xlabel('iteration',size=30)
    plt.xlabel('Delta t',size=30)
    plt.ylabel('Correlation',size=30)
    plt.xlim([0,2])
    plt.ylim([0,1])
    plt.show()
    
    ax = plt.subplot(111)
    pos = [0.2, 0.15, 0.7, 0.8]
    ax.set_position(pos)