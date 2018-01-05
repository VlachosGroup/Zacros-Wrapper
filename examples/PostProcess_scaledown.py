# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 19:12:15 2016

@author: mpnun
"""

import os
import sys
import numpy as np
import zacros_wrapper as zw

''' ------------ User input section ------------ '''
exe_file = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/stiff_input'
RunPath = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/ScaledownV3'
product_spec = 'B'                                  # product species
number_of_runs = 96
''' -------------------------------------------- '''

if __name__ == '__main__':                 # Need this line to make parallelization work

    # Count the iterations
    #n_folders = len([d for d in os.listdir(RunPath) if os.path.isdir(RunPath + d + '/')])
    n_folders = len(os.listdir(RunPath))

    print n_folders

    cum_batch = []
    for ind in range(1,n_folders+1):
        
        x = zw.Replicates()
        x.ParentFolder = os.path.join(RunPath, 'Iteration_' + str(ind))
        x.ReadMultipleRuns()

        if ind == 1:
            cum_batch = x
        else:
            cum_batch = zw.Replicates.time_sandwich(cum_batch, x)

    # Do numerical analysis
    cum_batch.ParentFolder = RunPath
    cum_batch.AverageRuns()
    cum_batch.ComputeStats(product_spec, window = [0.5, 1])
    cum_batch.WriteSA_output()
    cum_batch.PlotSensitivities()

    