# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 14:12:14 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling
from KMCrun import KMCrun
import numpy as np

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')
    
    # Set all directories
    exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
    KMC_source = 'C:/Users/mpnun/Desktop/NH3/'
    RunPath = 'C:/Users/mpnun/Desktop/NH3/Scaledown/'
    
    max_events = int(1e5)
    max_iterations = 15
    cutoff = 0.5
    ss_inc = 3.0
    write_summary = True
    n_samples = 100
    n_runs = 10
    n_procs = 4    
    Product = 'N2'    
    
    ''' Set up system '''
    y = KMCrun()
    y.data.Path = KMC_source
    y.data.ReadAllInput()
    y.exe_file = exe_file
    z = RateRescaling()
    z.batch.runtemplate = y
    z.scale_parent_fldr = RunPath
    
    ''' Run rescaling '''
    z.batch.ParentFolder = RunPath + 'Iteration_1/'
    z.batch.ReadMultipleRuns()
    z.batch.AverageRuns()
    delta_sdf = z.ProcessStepFreqs()         # compute change in scaledown factors based on simulation result
    
    # Change sampling
    z.batch.runtemplate.data.Conditions['MaxStep'] = 'inf'
    z.batch.runtemplate.data.Conditions['SimTime']['Max'] = z.batch.runAvg.data.Specnum['t'][-1]        
    z.batch.runtemplate.data.Report['procstat'] = ['time', z.batch.runAvg.data.Specnum['t'][-1] / n_samples]
    z.batch.runtemplate.data.Report['specnum'] = ['time', z.batch.runAvg.data.Specnum['t'][-1] / n_samples]
    
    is_steady_state = z.batch.runAvg.CheckSteadyState(Product)
    stiff = not np.max(np.abs(np.log10(delta_sdf))) < cutoff        # converged if changes to rate constants are small enough        
    
    # Update sampling
    if not is_steady_state:
        if stiff:
            z.batch.runtemplate.data.Conditions['SimTime']['Max'] = z.batch.runAvg.data.Specnum['t'][-1] / np.min(delta_sdf)
        else:
            z.batch.runtemplate.data.Conditions['SimTime']['Max'] = z.batch.runAvg.data.Specnum['t'][-1] * ss_inc        
        z.batch.runtemplate.data.Report['procstat'] = ['time', z.batch.runtemplate.data.Conditions['SimTime']['Max'] / n_samples]
        z.batch.runtemplate.data.Report['specnum'] = ['time', z.batch.runtemplate.data.Conditions['SimTime']['Max'] / n_samples]
    
    # Update the pre-exponential factors
    z.batch.runtemplate.AdjustPreExponentials(delta_sdf)
#    z.SDF_mat = np.vstack([z.SDF_mat, z.batch.runtemplate.data.scaledown_factors]) 
    
    
#    z.PerformScaledown(Product = 'B', n_runs = 10, n_procs = 4)
#    z.batch.runAvg.PlotElemStepFreqs()
#    z.PlotStiffnessReduction()
#    z.batch.runAvg.CheckSteadyState('B', show_graph = True)
#    print 'Final KMC time: ' + str(z.batch.runAvg.data.Specnum['t'][-1])