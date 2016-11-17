# -*- coding: utf-8 -*-
"""
Created on Fri Oct 07 19:12:15 2016

@author: mpnun
"""

import os
import sys

sys.path.insert(0, '../KMCsim')
from RateRescaling import RateRescaling

''' ------------ User input section ------------ '''
exe_file = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/zacros.exe'
KMC_source = 'C:/Users/mpnun/Desktop/Test/Source/'
RunPath = 'C:/Users/mpnun/Desktop/Test/RunDir/'
product_spec = 'B'                                  # product species
number_of_runs = 3
number_of_processors = 3
''' -------------------------------------------- '''

os.system('cls')

if __name__ == '__main__':                 # Need this line to make parallelization work

    z = RateRescaling()
    z.scale_parent_fldr = RunPath
    z.ReachSteadyState(product_spec, KMC_source, exe_file, max_iterations = 4, n_runs = number_of_runs, n_procs = number_of_processors)