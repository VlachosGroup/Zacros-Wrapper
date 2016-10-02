# -*- coding: utf-8 -*-
"""
Created on Sun Oct 02 13:22:25 2016

@author: mpnun
"""

import os
import sys

if __name__ == '__main__':                 # Need this line to make parallelization work

    sys.path.insert(0, '../KMCsim')
    from KMC_batch import KMC_batch
    from KMCrun import KMCrun
    
    os.system('cls')
    
    # Set all directories
    exe_path = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/'
    KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB_nonstiff/'
    RunPath = 'C:/Users/mpnun/Desktop/multirun_test/'
    
    # Read in data
    KMC_template = KMCrun()
    KMC_template.data.Path = KMC_source
    KMC_template.data.ReadAllInput()
    KMC_template.data.Path = RunPath
    KMC_template.exe_path = exe_path
    
    # Set up batch variables
    KMC_batch = KMC_batch()
    KMC_batch.runtemplate = KMC_template
    KMC_batch.ParentFolder = RunPath
    KMC_batch.n_runs = 4
    KMC_batch.Product = 'B'
    
    # Build folders and run jobs
    KMC_batch.BuildJobs()
    #KMC_batch.RunAllJobs(parallel = False)
    KMC_batch.RunAllJobs()
    
    # Analyze jobs
    KMC_batch.ReadMultipleRuns()
    KMC_batch.AverageRuns()
    print KMC_batch.runAvg.CheckSteadyState('B', show_graph = True)