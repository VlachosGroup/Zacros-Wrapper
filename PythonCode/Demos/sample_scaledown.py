# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 19:12:15 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling
from KMCrun import KMCrun

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')
    
    # Set all directories
    exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
    KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB/'
    RunPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/Tests/rescale_test/'
#    RunPath = 'C:/Users/mpnun/Desktop/'
    
    ''' Set up system '''
    y = KMCrun()
    y.data.Path = KMC_source
    y.data.ReadAllInput()
    y.exe_file = exe_file
    z = RateRescaling()
    z.batch.runtemplate = y
    z.scale_parent_fldr = RunPath
#    
    
    ''' Run rescaling '''
    z.PerformScaledown(Product = 'B', n_runs = 10, n_procs = 4)
    z.batch.runAvg.PlotElemStepFreqs()
#    z.ReadSummaryFile()
    
    ''' Analyze '''    
    z.PlotStiffnessReduction()
    z.PlotFinalTimes()
    z.batch.runAvg.CheckSteadyState('B', show_graph = True)