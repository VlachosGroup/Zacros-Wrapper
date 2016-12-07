# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import sys
from sys import getsizeof

sys.path.append('C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/PythonCode')
#sys.path.append('/home/vlachos/mpnunez/ZacrosWrapper')
import KMCsim as zw

################## User input ##################################

#zacros_exe = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
#KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/NonStiff/'
BatchPath = 'C:/Users/mpnun/Desktop/WGSsample/'
Product = 'CO2'
#n_runs = 10

################################################################

if __name__ == '__main__':                 # Need this line to make parallelization work

    # Batch of runs ----------------
    x = zw.Replicates()
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
    print getsizeof(x)

#    # Trajectory average
    x.AverageRuns()
#    
    x.runAvg.PlotSurfSpecVsTime()
#    x.runAvg.PlotGasSpecVsTime()
#    x.runAvg.PlotElemStepFreqs(window = [0.4, 1.0])
#    
#    x.runAvg.CalcRateTraj(Product)
#    x.runAvg.PlotRateVsTime()
#
##    
##     Rate and sensitivities
#    x.ComputeStats(Product, window = [0.4, 1])
#    x.PlotSensitivities()
#    x.WriteSA_output(BatchPath)