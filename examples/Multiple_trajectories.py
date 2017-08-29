# Read multiple trajectories and perform statistical analysis

import zacros_wrapper as zw

''' ------------ User input section ------------ '''
zacros_exe = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/Github/Zacros-Wrapper/input_files/AtoB/'
BatchPath = '/home/vlachos/mpnunez/test'
Product = 'B'
n_runs = 16
''' -------------------------------------------- '''

# Prepare template
x = zw.Replicates()
x.runtemplate = zw.kmc_traj()
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

# Perform sensitivity analysis
x.AverageRuns()
x.runAvg.PlotSurfSpecVsTime()
x.runAvg.PlotGasSpecVsTime()
x.runAvg.PlotElemStepFreqs(time_norm = True)