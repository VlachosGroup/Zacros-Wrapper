# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:26:19 2016

@author: mpnun
"""

from multiprocessing import Pool
import os
import numpy
import random
import sys

sys.path.insert(0, '../KMCsim')
from KMCrun import KMCrun
from AnalyzeData import AnalyzeData


os.system('cls')
exe_path = 'C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/Zacros_mod/'

def runKMC(KMC_inst):
    print KMC_inst.data.Path
    print '\n'
    KMC_inst.Run_sim(exe_path)
    return KMC_inst

if __name__ == '__main__':
    
    # Options     
#    run_dir = 'C:/Users/mpnun/Desktop/parallel_test/'
    run_dir = 'C:/Users/mpnun/Desktop/ss/' 
    
    n_runs = 50
    n_procs = 4
   
    # build KMC job from the input files
    KMC_template = KMCrun()
    KMC_template.data.Path = run_dir
    KMC_template.data.ReadAllInput()
   
    KMC_jobs = []
    for i in range(n_runs):
        new_job = KMCrun()
        new_job.data.Path = run_dir
        new_job.data.ReadAllInput()
        new_dir = run_dir + str(i+1) + '/'
        os.makedirs(new_dir)
        new_job.data.Path = new_dir
        
        new_job.data.Conditions['Seed'] = random.randint(10000, 99999)
        new_job.data.WriteAllInput()
        KMC_jobs.append(new_job)
        
    # Serial
#    for job in KMC_jobs:
#        runKMC(job)        
        
    # Parallel
    pool = Pool(processes=n_procs)
    finished_runs = pool.map(runKMC, KMC_jobs)
    
    # Analyze as a group of replicates
    run_group = AnalyzeData()
    run_group.ReadMultipleRuns(run_dir)
    run_group.AverageRuns()
    run_group.runAvg.PlotSurfSpecVsTime()
    run_group.runAvg.PlotElemStepFreqs()