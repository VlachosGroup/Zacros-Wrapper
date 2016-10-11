# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling
from KMCrun import KMCrun
from KMC_batch import KMC_batch

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')
    
    # Set all directories
    exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
    KMC_source = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/BigJobs/AtoB/'
    RunPath = 'C:/Users/mpnun/Desktop/rescale_test/'
    BatchPath = 'C:/Users/mpnun/Desktop/analyzethese/'


    # Test batch of runs ----------------
    x = KMC_batch()
    x.ParentFolder = BatchPath
    x.ReadMultipleRuns()
    x.AverageRuns()
    
    z = RateRescaling()
    z.batch = x
    z.batch.runtemplate = z.batch.runAvg
    
    delta_sdf = z.ProcessStepFreqs()
    print delta_sdf
    
#    z.batch.ParentFolder = 'C:/Users/mpnun/Desktop/COox/Scaledown/Iteration_2/'
#    z.batch.runtemplate.AdjustPreExponentials(delta_sdf)

# Need to finish this