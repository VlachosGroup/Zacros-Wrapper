# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from Replicates import Replicates

################## User input ##################################

BatchPath = 'C:/Users/mpnun/Documents/Local_research_files/ZacrosWrapper/AtoB_scaledown2/Iteration_4/'
#BatchPath = 'C:/Users/mpnun/Desktop/test_rep/'
Product = 'B'
n_cores = 3

################################################################

if __name__ == '__main__':                 # Need this line to make parallelization work

    os.system('cls')

    # Batch of runs ----------------
    x = Replicates()
    x.ParentFolder = BatchPath
    x.n_procs = n_cores  
    x.ReadMultipleRuns(parallel = True)

    # Trajectory average
    x.AverageRuns()
#    x.runAvg.PlotSurfSpecVsTime(save = False)
#    x.runAvg.PlotGasSpecVsTime(save = False)
#    x.runAvg.PlotElemStepFreqs(save = False)
    
    x.runAvg.CalcRateTraj(Product)
    x.runAvg.PlotRateVsTime()
    print x.runAvg.CheckSteadyState(Product, show_graph = True)
    print x.CheckAutocorrelation(Product)
    
#    x.runAvg.PlotPropsVsTime()
#    x.runAvg.PlotIntPropsVsTime()
    
#     Rate and sensitivities
#    x.ComputeStats(Product)
#    x.PlotSensitivities()
#    x.WriteSA_output(BatchPath)
    #x.WvarCheck()
#    
#    print [x.TOF, x.TOF_error]