# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 19:12:15 2016

@author: mpnun
"""

import os
import sys
import numpy as np

sys.path.append('/home/vlachos/mpnunez/ZacrosWrapper')
import KMCsim as zw

''' ------------ User input section ------------ '''
exe_file = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/stiff_input'
RunPath = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/ScaledownV3'
product_spec = 'B'                                  # product species
number_of_runs = 96
gas_stoich = np.array([-1, 1])
''' -------------------------------------------- '''

if __name__ == '__main__':                 # Need this line to make parallelization work
	
    z = zw.RateRescaling()
    z.scale_parent_fldr = RunPath
    z.ReachSteadyStateAndRescale(product_spec, gas_stoich, KMC_source, exe_file, n_runs = number_of_runs, platform = 'Squidward')