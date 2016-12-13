# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 17:58:26 2016

@author: mpnun
"""

# Analyze and make graphs

import os
import sys
import numpy as np

import matplotlib as mat
mat.use('Agg')
import matplotlib.pyplot as plt

sys.path.append('/home/1483/ZacrosWrapper')
import ZacrosWrapperPy as zw

################## User input ##################################

RunPath = '/home/1483/ZacrosWrapper/KMC_data/AtoB/Scaledown/'
Product = 'B'                                  # product species

################################################################

if __name__ == '__main__':                 # Need this line to make parallelization work

    # Count the iterations
    n_files = n_folders = 0
    for _, dirnames, filenames in os.walk(RunPath):
        n_files += len(filenames)
        n_folders += len(dirnames)

    iter_names = []
    time_vecs = []
    rate_vecs = []

    delt_list = []
    corr_list = []
    corr_err = []
    corr_expect_list = []
    for ind in range(1,n_folders+1):
        
        x = zw.Replicates()
        x.ParentFolder = RunPath + 'Iteration_' + str(ind) + '/'
        x.ReadMultipleRuns()
        
        # Trajectory average
        x.AverageRuns()
        
        x.runAvg.CalcRateTraj(Product)

        time_vecs.append(x.runAvg.Specnum['t'][1::])
        rate_vecs.append(x.runAvg.rate_traj)
        iter_names.append('Iteration_' + str(ind))

        corrs =  x.CheckAutocorrelation(Product)     
        
        delt = x.runAvg.Specnum['t'][-1]
        delt_list.append(delt)
        corr_list.append(corrs)
        corr_err.append(corrs[1])
        
 
    tau = 0.3186
    dtlin = np.linspace(0, delt_list[-1], 1000)
    corr_expected = []
    for dt in dtlin:
        corr_expected.append(np.exp(-dt / tau)) 
    
    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 12    
    
    plt.figure()
    plt.plot(dtlin, corr_expected)
    plt.errorbar(delt_list, corr_list, yerr = corr_err, markersize = 15, marker='o')
#    plt.plot(delt_list, corr_list, markersize = 15, marker='o')
    plt.xticks(size=24)
    plt.yticks(size=24)
#    plt.xlabel('iteration',size=30)
    plt.xlabel('Delta t',size=30)
    plt.ylabel('Correlation',size=30)
#    plt.xlim([0,7])
#    plt.ylim([0,1])
    plt.legend(['anal.', 'simulated'], loc=4, prop={'size':20}, frameon=False)
    ax = plt.subplot(111)
    pos = [0.2, 0.15, 0.7, 0.8]
    ax.set_position(pos)
    
    plt.savefig('corr_per_iter.png')
    plt.close()
    
    zw.Helper.PlotTrajectory(time_vecs, rate_vecs, series_labels = iter_names, xlab = 'time (s)', ylab = 'rate (1/s)', fname = 'rate_per_iter.png')