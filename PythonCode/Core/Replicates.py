# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 15:20:36 2016

@author: robieta
"""

import os
import numpy as np
import matplotlib as mat
mat.use('Agg')
import matplotlib.pyplot as plt
import copy

from KMC_Run import KMC_Run
from Helper import FileIO, Stats
import time

'''
Has data from muliple trajectories with the same input, but different 
random seeds and possibly different initial states
'''

class Replicates:
    
    def __init__(self):
             
        # General info
        self.ParentFolder = ''
        self.runtemplate = KMC_Run()                             # Use this to build replicate jobs
        self.runAvg = KMC_Run()            # Values are averages of all runs from runList
        self.n_runs = 0
        
        # Input data for different trajectories
        self.run_dirs = []
        self.rand_seeds = []
        
        # Output data taken accross all trajectories
        # Put here various multidimensional numpy arrays which will make it easier to do statistics
        self.t_vec = []                     # times at which data is recorded
        self.species_pops = []
        self.rxn_freqs = []
        self.History_final_snaps = []
        self.props = []
        self.Props_integ = []
        self.traj_derivs = []
        self.events_total = []
        self.CPU_total = []
        
        # Analysis
        self.Product = ''
        self.TOF = 0
        self.TOF_error = 0
        self.NSC_inst = []
        self.NSC_ci_inst = []
        self.NSC_erg = []
        self.NSC_ci_erg = []
    
    
    def BuildJobFiles(self, init_states = []):
        
        FileIO.ClearFolderContents(self.ParentFolder)    
        
        # List the directories and random seeds 
        self.run_dirs = []
        self.rand_seeds = []
        seed = 5000
        
        # Go through the directories and write input files for trajectories with
        # different random seeds and possibly different initial states
        for i in range(self.n_runs):
            
            # Set the random seed
            self.runtemplate.Conditions['Seed'] = seed
            self.rand_seeds.append(seed)
            seed = seed + 1
            
            # Set the path
            fldr = os.path.join(self.ParentFolder, str(i+1))
            self.runtemplate.Path = fldr
            self.run_dirs.append(fldr)
            
            # Set initial state
            if not init_states == []:
                self.runtemplate.StateInput['Type'] = 'history'
                self.runtemplate.StateInput['Struct'] = init_states[i]
                
            # Write input files in that folder
            if not os.path.exists(fldr):
                os.makedirs(fldr)
                
            self.runtemplate.WriteAllInput()
        
    def RunAllJobs_parallel_JobArray(self, max_cores = 100, server = 'Squidward'):
    
        with open(os.path.join(self.ParentFolder, 'dir_list.txt'), 'w') as txt:
            for fldr in self.run_dirs:
                txt.write(fldr + '\n')
        
        if server == 'Squidward':
            max_cores = np.min([max_cores, 96]) 
        elif server == 'Farber':
            max_cores = np.min([max_cores, 100]) 
        else:
            raise Exception('Unrecognized server for parallel runs')
        
        n_cores = np.min([max_cores, self.n_runs])           
        
        with open(os.path.join(self.ParentFolder, 'zacros_submit_JA.qs'), 'w') as txt:
            
            if server == 'Farber':             
            
                txt.write('#!/bin/bash\n')
                txt.write('#$ -cwd\n')
                txt.write('#$ -j y\n')
                txt.write('#$ -S /bin/bash\n')
                txt.write('#$ -l h_cpu=168:00:00\n')
                txt.write('#\n')
                txt.write('\n')
                txt.write('#$ -N zacros_JA 					#This is the name of the job array\n')
                txt.write('#$ -t 1-' + str(self.n_runs) + '  							#Assumes task IDs increment by 1; can also increment by another value\n')
                txt.write('#$ -tc ' + str(n_cores) + ' 							#This is the total number of tasks to run at any given moment\n')
                txt.write('#$ -pe threads 1 				#Change the last field to the number of processors desired per task\n')
                txt.write('#\n')
                txt.write('# Change the following to #$ and set the amount of memory you need\n')
                txt.write('# per-slot if you are getting out-of-memory errors using the\n')
                txt.write('# default:\n')
                txt.write('#$ -l m_mem_free=4G\n')
                txt.write('\n')
                txt.write('source /etc/profile.d/valet.sh\n')
                txt.write('\n')
                txt.write('# Use vpkg_require to setup the environment:\n')
                txt.write('vpkg_require intel/2016\n')
                txt.write('\n')
                txt.write('# Ensure that the OpenMP runtime knows how many processors to use;\n')
                txt.write('# Grid Engine automatically sets NSLOTS to the number of cores granted\n')
                txt.write('# to this job:\n')
                txt.write('export OMP_NUM_THREADS=$NSLOTS\n')
                txt.write('\n')
                txt.write('job_file=\'' + os.path.join(self.ParentFolder, 'dir_list.txt') + '\'\n')
                txt.write('#Change to the job directory\n')
                txt.write('job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")\n')
                txt.write('cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>\n')
                txt.write('                  #This could easily be modified to take a prefix; ask me how.\n')
                txt.write('\n')
                txt.write('# Now append whatever commands you use to run your OpenMP code:\n')
                txt.write(self.runtemplate.exe_file)
            
            else:       # Squidward
            
                txt.write('#!/bin/bash\n')
                txt.write('#$ -cwd\n')
                txt.write('#$ -j y\n')
                txt.write('#$ -S /bin/bash\n')
                txt.write('#\n')
                txt.write('\n')
                txt.write('#$ -N zacros_JA 					#This is the name of the job array\n')
                txt.write('#$ -t 1-' + str(self.n_runs) + '  							#Assumes task IDs increment by 1; can also increment by another value\n')
                txt.write('#$ -tc ' + str(n_cores) + ' 							#This is the total number of tasks to run at any given moment\n')
                txt.write('#$ -pe openmpi-smp 1 				#Change the last field to the number of processors desired per task\n')
                txt.write('\n')
                txt.write('job_file=\'' + os.path.join(self.ParentFolder, 'dir_list.txt') + '\'\n')
                txt.write('#Change to the job directory\n')
                txt.write('job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")\n')
                txt.write('cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>\n')
                txt.write('\n\n')
                txt.write(self.runtemplate.exe_file)
    
        # Call to system to submit the job array
        os.chdir(self.ParentFolder)
        os.system('qsub ' + os.path.join(self.ParentFolder, 'zacros_submit_JA.qs'))
    
        # Wait for jobs to be done
        all_jobs_done = False
        while not all_jobs_done:
            time.sleep(60)
            all_jobs_done = True
            for fldr in self.run_dirs:
                self.runtemplate.Path = fldr
                if not self.runtemplate.CheckComplete():
                    all_jobs_done = False
                    
        print 'Jobs in ' + self.ParentFolder + ' have finished'
    
    def RunAllJobs_serial(self):       # Serial version of running all jobs

        # loop over array of directories and execute the Zacros executable in that folder
        for fldr in self.run_dirs:
            self.runtemplate.Path = fldr
            self.runtemplate.Run_sim()
    
    def ReadMultipleRuns(self):     # Can take ~1 minutes to use this method
        
        print 'Reading all runs in ' + self.ParentFolder
        
        dummy_run = KMC_Run()       # Use this to transfer information
        
        # If directory list is empty, fill it with the directories in the current folder which have finished jobs
        if self.run_dirs == []:
        
            DirList = [d for d in os.listdir(self.ParentFolder) if os.path.isdir(os.path.join(self.ParentFolder, d))]      # List all folders in ParentFolder
            
            for direct in DirList:
            
                full_direct = os.path.join(self.ParentFolder, direct)
                dummy_run.Path = full_direct
                
                if dummy_run.CheckComplete():
                    self.run_dirs.append(full_direct)
                    
            self.n_runs = len(self.run_dirs)
        
        # Read first job to get some information
        dummy_run.Path = self.run_dirs[i]
        dummy_run.ReadAllOutput()
        
        # Create arrays for data
        # Use input data from runtemplate to properly size the arrays
        self.t_vec = np.array(dummy_run.Specnum['t'])
        self.n_timepoints = len(self.t_vec)
        
        n_specs = dummy_run.Species['n_gas'] + dummy_run.Species['n_surf']
        self.species_pops = np.zeros([self.n_runs, self.n_timepoints, n_specs])
        
        self.rxn_freqs = np.zeros([self.n_runs, self.n_timepoints, dummy_run.Reactions['nrxns']])
        self.History_final_snaps = []            # list of the final states
        self.props = np.zeros([self.n_runs, self.n_timepoints, dummy_run.Reactions['nrxns']])
        self.Props_integ = np.zeros([self.n_runs, self.n_timepoints, dummy_run.Reactions['nrxns']])
        self.traj_derivs = np.zeros([self.n_runs, self.n_timepoints, dummy_run.Reactions['nrxns']])
        self.events_total = np.zeros(self.n_runs)
        self.CPU_total = np.zeros(self.n_runs)
        
        for i in range(self.n_runs):
        
            # Switch to folder and read output files
            dummy_run.Path = self.run_dirs[i]
            dummy_run.ReadAllOutput()
            
            # Pass the data from this run into the larger data structures in 
            # so you will not have to store separate KMC_Run objects
            # The large arrays of data will be easier to process
            
            if not dummy_run.History == []:
                self.History_final_snaps.append(dummy_run.History[-1])
            self.events_total[i] = dummy_run.Performance['events_occurred']
            self.CPU_total[i] = dummy_run.Performance['CPU_time']
            
        print 'Total CPU time is ' + str(np.sum(self.CPU_total)) + ' seconds.'
    

    # Create a KMC run object with averaged species numbers, reaction firings, and propensities
    def AverageRuns(self):
        
        # Initialize run average with information from first run, then set data to zero
        self.runAvg = copy.deepcopy(self.runList[0])
        self.runAvg.Path = self.ParentFolder

        self.runAvg.Specnum['spec'] = np.zeros(self.runList[0].Specnum['spec'].shape)
        self.runAvg.Procstat['events'] = np.zeros(self.runList[0].Procstat['events'].shape)
        self.runAvg.Binary['propCounter'] = np.zeros(self.runList[0].Binary['propCounter'].shape)
        
        # Add data from each run
        for run in self.runList:
            self.runAvg.Specnum['spec'] = self.runAvg.Specnum['spec'] + run.Specnum['spec'].astype(float) / self.n_runs
            self.runAvg.Procstat['events'] = self.runAvg.Procstat['events'] + run.Procstat['events'].astype(float) / self.n_runs
            self.runAvg.Binary['propCounter'] = self.runAvg.Binary['propCounter'] + run.Binary['propCounter'] / self.n_runs

     
    def ComputeStats(self, product, window = [0.0, 1]):
        
        # Find indices for the beginning and end of the time window
        start_t = window[0] * self.runList[0].Specnum['t'][-1]
        end_t = window[1] * self.runList[0].Specnum['t'][-1]
        start_ind = self.runList[0].time_search(start_t)
        end_ind = self.runList[0].time_search(end_t)
        start_t = self.runList[0].Specnum['t'][start_ind]
        end_t = self.runList[0].Specnum['t'][end_t]        
        
        Tof_out = self.runAvg.ComputeTOF(product, win = window)
        tof_fracs_inst = Tof_out['TOF_fracs_inst']
        tof_fracs_inst = tof_fracs_inst[::2] + tof_fracs_inst[1::2]
        tof_fracs_erg = Tof_out['TOF_fracs_erg']
        tof_fracs_erg = tof_fracs_erg[::2] + tof_fracs_erg[1::2]
        
        Wdata = np.zeros([self.n_runs, self.runList[0].Reactions['nrxns']])      # number of runs x number of reactions
        rdata = np.zeros([self.n_runs, 2])     # 1st column: final rates, 2nd column: integral rates
        ind = 0
        for run in self.runList:
            Wdata[ind,:] = run.Binary['W_sen_anal'][end_ind,::2] - run.Binary['W_sen_anal'][start_ind,::2] + run.Binary['W_sen_anal'][end_ind,1::2] - run.Binary['W_sen_anal'][start_ind,1::2]
                               
            TOF_output = run.ComputeTOF(product, win = window)
            rdata[ind,0] = TOF_output['TOF_inst'] / Tof_out['TOF_inst']
            rdata[ind,1] = TOF_output['TOF_erg'] / Tof_out['TOF_erg']
            ind = ind + 1
        
#        # Rate data
#        TOF_stats = Stats.mean_ci(TOF_output)
#        self.TOF = TOF_stats[0]
#        self.TOF_error = TOF_stats[1]        
        
        # Sensitivity data
        A = np.hstack([Wdata, rdata])
        cov_out = Stats.cov_mat_ci( np.transpose(A) )    
        self.NSC_inst = cov_out['cov_mat'][:-2:,-2] + tof_fracs_inst
        self.NSC_ci_inst = cov_out['ci_mat'][:-2:,-2]
        self.NSC_erg = cov_out['cov_mat'][:-2:,-1] + tof_fracs_erg
        self.NSC_ci_erg = cov_out['ci_mat'][:-2:,-1]
        
        
                           
    def WvarCheck(self): 
        
        ''' Compute trajectory derivative variances vs. time '''        
        
        W_dims = self.runList[0].Binary['W_sen_anal'].shape
        n_timepoints = W_dims[0]
        n_rxns = W_dims[1]       
        
        Wvars = np.zeros([n_timepoints,n_rxns])
        for i in range(0,n_timepoints):
            for j in range(0,n_rxns):
                data_vec = np.zeros((self.n_runs))
                for k in range(0,self.n_runs):
                    data_vec[k] = self.runList[k].Binary['W_sen_anal'][i,j]
                Wvars[i,j] = np.var(data_vec)
        
        ''' Plot results '''
        
        FileIO.PlotOptions()
        plt.figure()
            
        labels = []
        for i in range (2*len(self.runList[0].Reactions['names'])):
            if np.max(np.abs( Wvars[:,i] )) > 0:
                plt.plot(self.runList[0].Specnum['t'], Wvars[:,i])
                labels.append(self.runList[0].Reactions['names'][i/2])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('Time (s)',size=24)
        plt.ylabel('var(W)',size=24)
        plt.legend(labels,loc=4,prop={'size':20},frameon=False)        
        plt.show()
        
    def PlotSensitivities(self): 
        
        FileIO.PlotOptions()
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        
        for i in range (self.runList[0].Reactions['nrxns']):
            cutoff = 0.05
            if self.NSC_inst[i] + self.NSC_ci_inst[i] > cutoff or self.NSC_inst[i] - self.NSC_ci_inst[i] < -cutoff:     
                plt.barh(ind-0.9, self.NSC_inst[i], width, color='r', xerr = self.NSC_ci_inst[i], ecolor='k')               
                ylabels.append(self.runList[0].Reactions['names'][i])              
                yvals.append(ind-0.6)                
                ind = ind - 1

        plt.plot([0, 0], [0, ind], color='k')
        plt.xlim([0,1])
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('NSC',size=24)
        plt.yticks(yvals, ylabels)
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        
        plt.savefig(os.path.join(self.ParentFolder, 'SA_inst_output.png'))
        plt.close()
        
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        
        for i in range (self.runList[0].Reactions['nrxns']):
            cutoff = 0.05
            if self.NSC_erg[i] + self.NSC_ci_erg[i] > cutoff or self.NSC_erg[i] - self.NSC_ci_erg[i] < -cutoff:     
                plt.barh(ind-0.9, self.NSC_erg[i], width, color='r', xerr = self.NSC_ci_erg[i], ecolor='k')               
                ylabels.append(self.runList[0].Reactions['names'][i])              
                yvals.append(ind-0.6)                
                ind = ind - 1

        plt.plot([0, 0], [0, ind], color='k')
        plt.xlim([0,1])
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('NSC',size=24)
        plt.yticks(yvals, ylabels)
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        
        plt.savefig(os.path.join(self.ParentFolder, 'SA_erg_output.png'))
        plt.close()
    
    def WriteSA_output(self):
        with open(os.path.join(self.ParentFolder, 'SA_output.txt'), 'w') as txt:
            txt.write('Normalized sensitivity coefficients \n\n')
            txt.write('Turnover frequency: \t' + '{0:.3E} \t'.format(self.TOF) + '+- {0:.3E} \t'.format(self.TOF_error) + '\n\n')               
            txt.write('Reaction name \t NSC \t NSC confidence \n')

            for rxn_ind in range(self.runList[0].Reactions['nrxns']):
                txt.write(self.runAvg.Reactions['names'][rxn_ind] + '\t' + '{0:.3f} +- \t'.format(self.NSC_inst[rxn_ind]) + '{0:.3f}\t'.format(self.NSC_ci_inst[rxn_ind]) + '{0:.3f} +- '.format(self.NSC_erg[rxn_ind]) + '{0:.3f}'.format(self.NSC_ci_erg[rxn_ind]) + '\n')
                
    def FD_SA(self, rxn_inds = [1], pert_frac = 0.05, n_runs = 20, setup = True, exec_run = True, analyze_bool = True):
        
        # Create objects for perturbed systems
        plus = copy.deepcopy(self)
        minus = copy.deepcopy(self)
        FD_list = [plus, minus]

        for FD in FD_list:
            FD.n_runs = n_runs
        
        # Adjust pre-exponential factors in each
        adjust_plus = np.ones(self.runtemplate.Reactions['nrxns'])
        adjust_minus = np.ones(self.runtemplate.Reactions['nrxns'])
        for rxn_ind in rxn_inds:
            adjust_plus[rxn_ind-1] = 1 + pert_frac
            adjust_minus[rxn_ind-1] = 1 / (1 + pert_frac)
        plus.runtemplate.AdjustPreExponentials(adjust_plus)
        minus.runtemplate.AdjustPreExponentials(adjust_minus)
        
        # Set subfolder for perturbed runs
        plus.ParentFolder = self.ParentFolder + 'plus'
        minus.ParentFolder = self.ParentFolder + 'minus'
        
        if setup:
        
            for FD in FD_list:
                for rxn_ind in rxn_inds:
                    FD.ParentFolder = FD.ParentFolder + '_' + str(rxn_ind)
                FD.ParentFolder = FD.ParentFolder + '/'
            
                # Build folders for runs
                if not os.path.exists(FD.ParentFolder):
                        os.makedirs(FD.ParentFolder)
                FD.BuildJobs()
        
        if exec_run:
            for FD in FD_list:
                FD.RunAllJobs()
        
        ''' Analyze results '''
        if analyze_bool:
            for FD in FD_list:
                FD.ReadMultipleRuns()
                FD.AverageRuns()
                
            plus_stats = []
            for run in plus.runList:
                plus_stats.append(run.ComputeTOF(self.Product)['TOF'])
                
            minus_stats = []
            for run in minus.runList:
                minus_stats.append(run.ComputeTOF(self.Product)['TOF'])
    
            all_TOFs = plus_stats + minus_stats
            TOF_mean = np.mean(all_TOFs)
            diff_stats = Stats.diff_ci(plus_stats, minus_stats)   
            
            NSC = diff_stats[0] / TOF_mean / (2 * pert_frac)
            NSC_ci = diff_stats[1] / TOF_mean / (2 * pert_frac)
            
            return [NSC, NSC_ci]
            
    def CheckAutocorrelation(self, Product, limits = [0.0, 1]):
        
        data1 = []
        data2 = []
        for run in self.runList:
            run.CalcRateTraj(Product)
            
#            ind1 = run.fraction_search(limits[0])
#            ind2 = run.fraction_search(limits[1])
            
            data1.append(run.rate_traj[0])
            data2.append(run.rate_traj[-1])
        
        var1 = np.var(data1)
        var2 = np.var(data2)
        if var1 == 0 or var2 == 0:
            return [1.0, 0.0]
        else:
            return Stats.cov_ci(data1,data2) / np.sqrt( np.var(data1) * np.var(data2) )
            
    @staticmethod
    def time_sandwich(batch1, batch2):
        
        sand = copy.deepcopy(batch2)        
        
        if batch1.n_runs == 1:
            return sand
        
        for run_ind in range(batch1.n_runs):
            sand.runList[run_ind] = KMC_Run.time_sandwich(batch1.runList[run_ind], sand.runList[run_ind])
            
        return sand