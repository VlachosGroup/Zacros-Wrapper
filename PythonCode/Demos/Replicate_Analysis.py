# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import sys

from mpi4py import MPI

sys.path.insert(0, '../KMCsim')
from Replicates import Replicates

################## User input ##################################

zacros_exe = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/NonStiff/'
BatchPath = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/test_parallel/'
Product = 'B'
n_runs = 10

################################################################

if __name__ == '__main__':                 # Need this line to make parallelization work

    # Batch of runs ----------------
    x = Replicates()
    x.ParentFolder = BatchPath
#    x.runtemplate.Path = KMC_source
#    x.runtemplate.exe_file = zacros_exe
#    x.runtemplate.ReadAllInput()    
    
    # Build and run
#    x.n_runs = n_runs
#    x.BuildJobsFromTemplate()
#    x.BuildJobFiles()
#    x.RunAllJobs()
    
    x.ReadMultipleRuns()

    # Trajectory average
    x.AverageRuns()
    
    COMM = MPI.COMM_WORLD
    if COMM.rank == 0:
        x.runAvg.PlotSurfSpecVsTime()
        x.runAvg.PlotGasSpecVsTime()
        x.runAvg.PlotElemStepFreqs()
    
#    x.runAvg.CalcRateTraj(Product)
#    x.runAvg.PlotRateVsTime()
#    print x.runAvg.CheckSteadyState(Product, show_graph = True)
#    print x.CheckAutocorrelation(Product)
#    
#    x.runAvg.PlotPropsVsTime()
#    x.runAvg.PlotIntPropsVsTime()
#    
#     Rate and sensitivities
#    x.ComputeStats(Product)
#    x.PlotSensitivities()
#    x.WriteSA_output(BatchPath)
#    x.WvarCheck()
    
#    print [x.TOF, x.TOF_error]