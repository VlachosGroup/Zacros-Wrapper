# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 19:12:15 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling
from KMC_Run import KMC_Run

''' ------------ User input section ------------ '''
exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB/'
RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_ZW/Scaledown/'
product_spec = 'B'                                  # product species
number_of_runs = 10
number_of_processors = 4
''' -------------------------------------------- '''

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')		# Clear the screen
	
    ''' Set up system '''
    y = KMC_Run()
    y.Path = KMC_source
    y.ReadAllInput()
    y.exe_file = exe_file
    z = RateRescaling()
    z.batch.runtemplate = y
    z.scale_parent_fldr = RunPath
#    
    
    ''' Run rescaling '''
    z.PerformScaledown(Product = product_spec, n_runs = number_of_runs, n_procs = number_of_processors)
    z.batch.runAvg.PlotElemStepFreqs()
#    z.ReadSummaryFile()            # Use this line if it has already been run
    
    ''' Analyze '''    
    z.PlotStiffnessReduction()
    z.PlotFinalTimes()
    z.batch.runAvg.CheckSteadyState(product_spec, show_graph = True)