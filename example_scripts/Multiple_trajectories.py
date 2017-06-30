# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 13:48:34 2016

@author: mpnun
"""

import sys
sys.path.append('/home/vlachos/mpnunez/Github/Zacros-Wrapper-expt/Zacros-Wrapper/PythonCode')
import Core as zw

################## User input ##################################

zacros_exe = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/KMC_data/AtoB/NonStiff/'
BatchPath = '/home/vlachos/mpnunez/ZacrosWrapper/KMC_data/AtoB/test'
#Product = 'B'
n_runs = 16

################################################################

# Prepare template
x = zw.Replicates()
x.runtemplate.Path = KMC_source
x.runtemplate.exe_file = zacros_exe
x.runtemplate.ReadAllInput()
x.ParentFolder = BatchPath

# Build files and run
x.n_runs = n_runs
x.BuildJobFiles()
x.RunAllJobs_parallel_JobArray()

# Read results
x.ReadMultipleRuns()

# Perform analysis