# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:58:26 2016

@author: mpnun
"""

# Analyze and make graphs

import os
import sys

import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np
import copy

sys.path.insert(0, '../KMCsim')
from Replicates import Replicates
from KMC_Run import KMC_Run

################## User input ##################################

RunPath = 'C:/Users/mpnun/Desktop/WGS_ss/'
Product = 'CO2'
n_iterations = 6

################################################################

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')

    cum_batch = Replicates()

    delt_list = []
    corr_list = []
    corr_err = []
    corr_expect_list = []
    for ind in range(2,n_iterations+1):
        
        x = Replicates()
        x.ParentFolder = RunPath + 'Iteration_' + str(ind) + '/'
        x.ReadMultipleRuns()
        
        # Add data to running list
        if ind == 1:
            pass            # Do not use data from first run because it is on event rather than time
        elif ind == 2:
            cum_batch = copy.deepcopy(x)
        elif ind > 2:
            for run_ind in range(x.n_runs):
                cum_batch.runList[run_ind] = KMC_Run.time_sandwich(cum_batch.runList[run_ind], x.runList[run_ind])        
        
        # Trajectory average
        cum_batch.AverageRuns()
        
        cum_batch.runAvg.CalcRateTraj(Product)
#        print x.runAvg.CheckSteadyState(Product)
        corrs =  cum_batch.CheckAutocorrelation(Product)     
        
        delt = cum_batch.runAvg.Specnum['t'][-1]/2
        delt_list.append(delt)
        corr_list.append(corrs)
#        corr_err.append(corrs[1])
        
 
#    tau = 0.3186
#    dtlin = np.linspace(0, delt_list[-1], 1000)
#    corr_expected = []
#    for dt in dtlin:
#        corr_expected.append(np.exp(-dt / tau)) 
    
    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 12    
    
    plt.figure()
#    plt.errorbar(range(2,8), corr_list, yerr = corr_err, markersize = 15, marker='o')
#    plt.plot(dtlin, corr_expected)
#    plt.errorbar(delt_list, corr_list, yerr = corr_err, markersize = 15, marker='o')
    plt.plot(delt_list, corr_list, markersize = 15, marker='o')
    plt.xticks(size=24)
    plt.yticks(size=24)
#    plt.xlabel('iteration',size=30)
    plt.xlabel('Delta t',size=30)
    plt.ylabel('Correlation',size=30)
#    plt.xlim([0,7])
#    plt.ylim([0,1])
    plt.show()
    
    ax = plt.subplot(111)
    pos = [0.2, 0.15, 0.7, 0.8]
    ax.set_position(pos)