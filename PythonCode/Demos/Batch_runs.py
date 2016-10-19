# -*- coding: utf-8 -*-
"""
Created on Sun Oct 02 13:22:25 2016

@author: mpnun
"""

import os
import sys
import numpy as np

if __name__ == '__main__':                 # Need this line to make parallelization work

    sys.path.insert(0, '../KMCsim')
    from KMC_batch import KMC_batch
    from KMCrun import KMCrun
    
    os.system('cls')
    
    # Set all directories
    exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
    KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_ZW/LRSA/1/'
    RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_ZW/FDSA/'
    
    # Read in data
    KMC_template = KMCrun()
    KMC_template.data.Path = KMC_source
    KMC_template.data.ReadAllInput()
    KMC_template.data.Path = RunPath
    KMC_template.exe_file = exe_file
    KMC_template.data.scaledown_factors = np.ones(KMC_template.data.Reactions['nrxns'])
    
    # Set up batch variables
    KMC_batch = KMC_batch()
    KMC_batch.runtemplate = KMC_template
    KMC_batch.ParentFolder = RunPath
    KMC_batch.n_runs = 20
    KMC_batch.n_procs = 3
    KMC_batch.Product = 'B'
    
    NSC_1 = KMC_batch.FD_SA(rxn_inds = [1])
    NSC_2 = KMC_batch.FD_SA(rxn_inds = [2])
    NSC_3 = KMC_batch.FD_SA(rxn_inds = [3])
    
    print [NSC_1[0], NSC_2[0], NSC_3[0]]
    print [NSC_1[1], NSC_2[1], NSC_3[1]]
    
#    # Build folders and run jobs
#    KMC_batch.BuildJobs()
#    KMC_batch.RunAllJobs()
#    
#    # Analyze jobs
#    KMC_batch.ReadMultipleRuns()
#    KMC_batch.AverageRuns()
#    print KMC_batch.runAvg.CheckSteadyState('B', show_graph = True)