# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 19:12:15 2016

@author: mpnun
"""

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mat

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling
from KMC_Run import KMC_Run

''' ------------ User input section ------------ '''
KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB/'
RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_scaledown2/'
product_spec = 'B'                                  # product species
number_of_runs = 10
number_of_processors = 4
''' -------------------------------------------- '''

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')		# Clear the screen
	
    n_iters = len([d for d in os.listdir(RunPath) if os.path.isdir(RunPath + d + '/')])
    iterations = np.zeros([n_iters,1])
    performance_data = np.zeros([n_iters,3]) 
    
    for i in range(n_iters):
        iterations[i,0] = i+1
        y = KMC_Run()
        y.Path = RunPath + 'Iteration_' + str(i+1) + '/1/'
        y.ReadSimIn()
        y.ReadGeneral()
        performance_data[i,:] = [y.Performance['t_final'], y.Performance['events_occurred'], y.Performance['CPU_time']]
        
    t_final_converged = performance_data[-1,0]
    CPU_required = np.zeros([n_iters,1])
    for i in range(n_iters):
        CPU_required[i,0] = t_final_converged / performance_data[i,0] * performance_data[i,2]
    
    # Plotting
    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 16
    
    plt.figure()
    
    plt.plot(iterations[:,0], CPU_required[:,0] / 3600, 'o-', markersize = 15)
    
    plt.xticks(size=24)
    plt.yticks(size=24)
    plt.xlabel('iterations',size=30)
    plt.ylabel('CPU requirement (hours)',size=30)
    plt.xlim([0,n_iters])
    plt.show()
    
    plt.yscale('log')
    ax = plt.subplot(111)
    pos = [0.2, 0.15, 0.7, 0.8]
    ax.set_position(pos)