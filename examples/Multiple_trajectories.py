# Read multiple trajectories and perform statistical analysis

import zacros_wrapper as zw

''' ------------ User input section ------------ '''
zacros_exe = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/KMC_data/AtoB/NonStiff/'
BatchPath = '/home/vlachos/mpnunez/ZacrosWrapper/KMC_data/AtoB/test'
#Product = 'B'
n_runs = 16
''' -------------------------------------------- '''

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

# Perform sensitivity analysis
# ...