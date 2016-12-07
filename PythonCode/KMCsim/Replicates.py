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
from Stats import Stats
from Helper import Helper
import time

class Replicates:
    
    def __init__(self):
             
        # General info
        self.ParentFolder = ''
        self.runList = []              # List of KMC_Run objects from which data is averaged
        self.runtemplate = KMC_Run()                             # Use this to build replicate jobs
        self.runAvg = KMC_Run()            # Values are averages of all runs from runList
        self.n_runs = 0
        
        # Analysis
        self.Product = ''
        self.TOF = 0
        self.TOF_error = 0
        self.NSC_inst = []
        self.NSC_ci_inst = []
        self.NSC_erg = []
        self.NSC_ci_erg = []
    
    def BuildJobsFromTemplate(self):
        
        # Build list of KMC_Run objects
        self.runList = []
        for i in range(self.n_runs):
            new_run = copy.deepcopy(self.runtemplate)
            new_run.Conditions['Seed'] = self.runtemplate.Conditions['Seed'] + i
            new_run.Path = self.ParentFolder + str(i+1) + '/'
            self.runList.append(new_run)
            
    def BuildJobFiles(self, write_dir_list = True, max_cores = 100):
        
        Helper.ClearFolderContents(self.ParentFolder)    
        
        # Build folders and input files for each job
        for run in self.runList:
            if not os.path.exists(run.Path):
                os.makedirs(run.Path)
            run.WriteAllInput()
            
        if write_dir_list:
            with open(self.ParentFolder + 'dir_list.txt', 'w') as txt:
                for job in self.runList:
                    txt.write(job.Path + '\n')
            
            n_cores = np.min([max_cores, self.n_runs])
            with open(self.ParentFolder + 'zacros_submit_JA.qs', 'w')as txt:
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
                txt.write('job_file=\'' + self.ParentFolder + 'dir_list.txt\'\n')
                txt.write('#Change to the job directory\n')
                txt.write('job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")\n')
                txt.write('cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>\n')
                txt.write('                  #This could easily be modified to take a prefix; ask me how.\n')
                txt.write('\n')
                txt.write('# Now append whatever commands you use to run your OpenMP code:\n')
                txt.write(self.runtemplate.exe_file)
    
    def SubmitJobArray(self):
        os.chdir(self.ParentFolder)
        os.system('qsub ' + self.ParentFolder + 'zacros_submit_JA.qs')
    
    def WaitForJobs(self):
        
        all_jobs_done = False
        while not all_jobs_done:
            time.sleep(60)
            all_jobs_done = True
            for job in self.runList:
                if not job.CheckComplete():
                    all_jobs_done = False
    
    def RunAllJobs(self):

        for job in self.runList:
            job.Run_sim()
    
    def ReadMultipleRuns(self):

        self.runList = []
        DirList = [d for d in os.listdir(self.ParentFolder) if os.path.isdir(self.ParentFolder + d + '/')]      # List all folders in ParentFolder
        for direct in DirList:
            run = KMC_Run()
            run.Path =  self.ParentFolder + direct + '/'
            if run.CheckComplete():
                self.runList.append(run)
        self.n_runs = len(self.runList)

        for job in self.runList:
            job.ReadAllOutput()
    
    @staticmethod
    def ReadPerformance(path):
        
        t_final_cum = 0
        events_occurred_cum = 0
        CPU_time_cum = 0        
        n_runs = 0        
        
        DirList = [d for d in os.listdir(path) if os.path.isdir(path + d + '/')]      # List all folders in ParentFolder
        for direct in DirList:
            run = KMC_Run()
            run.Path =  path + direct + '/'
            if run.CheckComplete():
                run.ReadAllInput()
                run.ReadGeneral()
                t_final_cum += run.Performance['t_final']
                events_occurred_cum += run.Performance['events_occurred']
                CPU_time_cum += run.Performance['CPU_time']
                n_runs += 1
                
        with open(path + 'Performance_summary.txt', 'w') as txt:   
            txt.write('----- Performance totals -----\n')
            txt.write('number of runs: ' + str(n_runs) + '\n' )      # number of runs
            txt.write('KMC time: {0:.3E} \n'.format(t_final_cum))      # seconds
            txt.write('events: ' + str(events_occurred_cum) + '\n' )                                  # events
            txt.write('CPU time: {0:.3E} \n'.format(CPU_time_cum))                       # seconds
            
        return {'t_final_cum': t_final_cum, 'events_occurred_cum': events_occurred_cum, 'CPU_time_cum': CPU_time_cum, 'n_runs': n_runs}

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

     
    def ComputeStats(self, product, window = [0.5, 1]):
        
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
                               
            TOF_output = run.ComputeTOF(product)
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
        
        Helper.PlotOptions()
        plt.figure()
            
        labels = []
        for i in range (2*len(self.runList[0].Reactions['names'])):
            if np.max(np.abs( Wvars[:,i] )) > 0:
                plt.plot(self.runList[0].Specnum['t'], Wvars[:,i])
                labels.append(self.runList[0].Reactions['names'][i/2])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('var(W)',size=24)
        plt.legend(labels,loc=4,prop={'size':20},frameon=False)        
        plt.show()
        
    def PlotSensitivities(self): 
        
        Helper.PlotOptions()
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
        
        plt.savefig(self.ParentFolder + 'SA_inst_output.png')
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
        
        plt.savefig(self.ParentFolder + 'SA_erg_output.png')
        plt.close()
    
    def WriteSA_output(self):
        with open(self.ParentFolder + 'SA_output.txt', 'w') as txt:
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
            
    def CheckAutocorrelation(self, Product, limits = [0.5, 1]):
        
        data1 = []
        data2 = []
        for run in self.runList:
            run.CalcRateTraj(Product)
            
            ind1 = run.time_search(run.Specnum['t'][-1] * limits[0])
            ind2 = run.time_search(run.Specnum['t'][-1] * limits[1])
            
            data1.append(run.rate_traj[ind1-1])
            data2.append(run.rate_traj[ind2-1])
        
        return Stats.cov_ci(data1,data2) / np.var(data2)
        #return Stats.cov_calc(data1,data2) / np.var(data2)