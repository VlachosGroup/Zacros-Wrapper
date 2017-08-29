# Read multiple trajectories and perform statistical analysis

import zacros_wrapper as zw
from mpi4py import MPI

''' ------------ User input section ------------ '''
zacros_exe = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/Github/Zacros-Wrapper/input_files/AtoB/'
BatchPath = '/home/vlachos/mpnunez/test'
Product = 'B'
n_runs = 16
''' -------------------------------------------- '''

if __name__ == '__main__':

    COMM = MPI.COMM_WORLD
    COMM.Barrier()

    # Prepare template
    x = zw.Replicates()
    x.runtemplate = zw.kmc_traj()
    x.runtemplate.Path = KMC_source
    x.runtemplate.exe_file = zacros_exe
    x.runtemplate.ReadAllInput()
    x.ParentFolder = BatchPath
    
    # Build files and run
    x.n_runs = n_runs
    if COMM.rank == 0:
        x.BuildJobFiles()
    
    '''
    Run all trajectories using MPI parallelization
    '''

    # Collect whatever has to be done in a list. Here we'll just collect a list of
    # numbers. Only the first rank has to do this.
    if COMM.rank == 0:
        jobs = x.run_dirs
        jobs = [jobs[_i::COMM.size] for _i in range(COMM.size)]             # Split into however many cores are available.
    else:
        jobs = None
    
    jobs = COMM.scatter(jobs, root=0)           # Scatter jobs across cores.
    
    # Now each rank just does its jobs and collects everything in a results list.
    # Make sure to not use super big objects in there as they will be pickled to be
    # exchanged over MPI.
    for job in jobs:
        x.runtemplate.Path = job
        x.runtemplate.Run_sim() 
    
    '''
    Read and analyze results
    '''
    if COMM.rank == 0:
        x.ReadMultipleRuns()
        x.AverageRuns()
        x.runAvg.PlotSurfSpecVsTime()
        x.runAvg.PlotGasSpecVsTime()
        x.runAvg.PlotElemStepFreqs(time_norm = True)