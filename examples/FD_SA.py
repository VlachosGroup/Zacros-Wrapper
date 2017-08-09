# -*- coding: utf-8 -*-
"""
Created on Sun Oct 02 13:22:25 2016

@author: mpnun
"""

import os
import sys
import numpy as np


from Replicates import Replicates
from KMC_Run import KMC_Run

''' ------------ User input section ------------ '''
exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_ZW/LRSA/1/'
RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_ZW/FDSA/'
number_of_runs = 20
number_of_processors = 3
product_species = 'B'
''' -------------------------------------------- '''

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')
        
    # Read in data
    KMC_template = KMC_Run()
    KMC_template.Path = KMC_source
    KMC_template.ReadAllInput()
    KMC_template.Path = RunPath
    KMC_template.exe_file = exe_file
    KMC_template.scaledown_factors = np.ones(KMC_template.Reactions['nrxns'])
    
    # Set up batch variables
    Replicates = Replicates()
    Replicates.runtemplate = KMC_template
    Replicates.ParentFolder = RunPath
    Replicates.n_runs = number_of_runs
    Replicates.n_procs = number_of_processors
    Replicates.Product = product_species
    
    NSC_1 = Replicates.FD_SA(rxn_inds = [1])
    NSC_2 = Replicates.FD_SA(rxn_inds = [2])
    NSC_3 = Replicates.FD_SA(rxn_inds = [3])
    
    print [NSC_1[0], NSC_2[0], NSC_3[0]]
    print [NSC_1[1], NSC_2[1], NSC_3[1]]
    
#    # Build folders and run jobs
#    Replicates.BuildJobs()
#    Replicates.RunAllJobs()
#    
#    # Analyze jobs
#    Replicates.ReadMultipleRuns()
#    Replicates.AverageRuns()
#    print Replicates.runAvg.CheckSteadyState('B', show_graph = True)