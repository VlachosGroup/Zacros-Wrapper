# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import sys
import numpy as np

sys.path.append('C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/PythonCode')
#sys.path.append('/home/vlachos/mpnunez/ZacrosWrapper')
import KMCsim as zw

if __name__ == '__main__':                 # Need this line to make parallelization work

    # Batch of runs ----------------
    b1 = zw.Replicates()
    b1.ParentFolder = 'C:\Users\mpnun\Documents\Local_research_files\ZacrosWrapper\iter_example\1'
    b1.ReadMultipleRuns()
    
    b2 = zw.Replicates()
    b2.ParentFolder = 'C:\Users\mpnun\Documents\Local_research_files\ZacrosWrapper\iter_example\2'
    b2.ReadMultipleRuns()
    
    b3 = zw.Replicates()
    b3.ParentFolder = 'C:\Users\mpnun\Documents\Local_research_files\ZacrosWrapper\iter_example\3'
    b3.ReadMultipleRuns()
    
    b2 = zw.Replicates.time_sandwich(b1, b2)
    b3 = zw.Replicates.time_sandwich(b2, b3)

    b3.AverageRuns()
#    
    b3.runAvg.PlotSurfSpecVsTime()
    b3.runAvg.PlotGasSpecVsTime()
    
    
    b3.runAvg.gas_stoich = np.array([-1, 0, 1, 1, -1])        # stoichiometry of gas-phase reaction
    b3.runAvg.calc_net_rxn()
    b3.runAvg.PlotNetGasRxnVsTime()
    print b3.runAvg.CheckNetRxnConvergence()
    
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