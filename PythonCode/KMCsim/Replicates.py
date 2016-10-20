# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 15:20:36 2016

@author: robieta
"""

import os, shutil
import pickle
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import copy

from KMC_Run import KMC_Run
import Helper as ut
from Stats import Stats

# Separate function for use in parallelization
def runKMC(kmc_rep):
    kmc_rep.Run_sim()

class Replicates:
    
    def __init__(self):
             
        # General info
        self.ParentFolder                     = ''
        self.runList                          = []
        self.runtemplate = KMC_Run()                             # Use this to build replicate jobs
        self.runAvg                           = KMC_Run()            # Values are averages of all runs
        self.n_runs = 0               
        self.n_procs = 4        
        
        # Analysis
        self.Product                          = ''
        self.TOF                              = 0
        self.TOF_error                        = 0
        self.NSC                              = []
        self.NSC_ci                           = []
    
    def BuildJobs(self):
        
        # Build list of KMC_Run objects
        self.runList = []
        for i in range(self.n_runs):
            new_run = copy.deepcopy(self.runtemplate)
            new_run.data.Conditions['Seed'] = self.runtemplate.data.Conditions['Seed'] + i
            new_run.data.Path = self.ParentFolder + str(i+1) + '/'
            self.runList.append(new_run)
        
        # Delete all files and folders in the run directory
        for the_file in os.listdir(self.ParentFolder):
            file_path = os.path.join(self.ParentFolder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)
                
        # Build folders and input files for each job
        for run in self.runList:
            if not os.path.exists(run.data.Path):
                os.makedirs(run.data.Path)
            run.data.WriteAllInput()
    
    def RunAllJobs(self, parallel = True):

        if parallel:
            pool = Pool(processes = self.n_procs)
            pool.map(runKMC, self.runList)
            pool.close()
        else:
            for run in self.runList:
                run.Run_sim()
    
    def ReadMultipleRuns(self):     # Could simplify this
        
#        print 'Reading output from ' + self.ParentFolder
        summary_fname = 'BatchSummary.p'
        if os.path.isfile(self.ParentFolder + summary_fname):                              # Runs have already been read
            print 'Existing pickle file found.'
            self.runList = pickle.load(open( self.ParentFolder + summary_fname, "rb" ))
        else:                                                               # Go through each directory and read the data           
            DirList = ut.Helper().GetDir(self.ParentFolder)
            nDir = len(DirList)
            self.runList = ['' for i in range(nDir)]
            for i in range(nDir):
#                print 'Reading run # ' + str(i+1) + ' / ' + str(nDir)
                RunPath = self.ParentFolder + DirList[i] + '/'
                self.runList[i]  = KMC_Run()
                self.runList[i].data.Path = RunPath
                self.runList[i].data.ReadAllOutput()                               # Read input and output files
            pickle.dump( self.runList, open( self.ParentFolder + summary_fname, "wb" ) )  
        self.n_runs = len(self.runList)

    # Create a KMC run object with averaged species numbers, reaction firings, and propensities
    def AverageRuns(self):
        
        # Initialize run average with information from first run, then set data to zero
        self.runAvg = KMC_Run()      
        self.runAvg.data = self.runList[0].data
        self.runAvg.data.Specnum['t'] = self.runList[0].data.Specnum['t']
             
        self.runAvg.data.Specnum['spec'] = self.runList[0].data.Specnum['spec'] - self.runList[0].data.Specnum['spec']         
        self.runAvg.data.Procstat['events'] = self.runList[0].data.Procstat['events'] - self.runList[0].data.Procstat['events']
#        self.runAvg.data.Binary['cluster'] = self.runList[0].data.Binary['cluster'] - self.runList[0].data.Binary['cluster']
#        self.runAvg.data.Binary['prop'] = self.runList[0].data.Binary['prop'] - self.runList[0].data.Binary['prop']
        self.runAvg.data.Binary['propCounter'] = self.runList[0].data.Binary['propCounter'] - self.runList[0].data.Binary['propCounter']        
        
        self.runAvg.data.Specnum['spec'] = self.runAvg.data.Specnum['spec'].astype(float)     
        self.runAvg.data.Procstat['events'] = self.runAvg.data.Procstat['events'].astype(float) 
#        self.runAvg.data.Binary['cluster'] = self.runAvg.data.Binary['cluster'].astype(float)         
        
        # Add data from each run
        for run in self.runList:
            self.runAvg.data.Specnum['spec'] = self.runAvg.data.Specnum['spec'] + run.data.Specnum['spec'].astype(float) / self.n_runs     
            self.runAvg.data.Procstat['events'] = self.runAvg.data.Procstat['events'] + run.data.Procstat['events'].astype(float) / self.n_runs
#            self.runAvg.data.Binary['cluster'] = self.runAvg.data.Binary['cluster'] + run.data.Binary['cluster'].astype(float) / self.n_runs
#            self.ru1nAvg.data.Binary['prop'] = self.runAvg.data.Binary['prop'] + run.data.Binary['prop'] / self.n_runs
            self.runAvg.data.Binary['propCounter'] = self.runAvg.data.Binary['propCounter'] + run.data.Binary['propCounter'] / self.n_runs

     
    def ComputeStats(self, product, SA = True):
        
        Tof_out = self.runAvg.ComputeTOF(product)
        tof_fracs = Tof_out['TOF_fracs']          
        
        TOF_vec = []
        for run in self.runList:
            TOF_vec.append(run.ComputeTOF(product)['TOF'])
        
        self.TOF = Stats.mean_ci(TOF_vec)[0]
        self.TOF_error = Stats.mean_ci(TOF_vec)[1]   
        
        if not SA:
            return        
        
        Wdata = np.zeros([self.n_runs, 2*self.runList[0].data.Reactions['nrxns']])      # number of runs x number of reactions
        TOFdata = np.zeros(self.n_runs)
        ind = 0
        for run in self.runList:
            Wdata[ind,:] = run.data.Binary['W_sen_anal'][-1,:]
#            Wdata[ind,:] = run.data.Procstat['events'][-1,:] - run.data.Binary['propCounter'][-1,:]
                               
            TOF_output = run.ComputeTOF(product)
            TOFdata[ind] = TOF_output['TOF']
            ind = ind + 1
        
        self.NSC = np.zeros(self.runList[0].data.Reactions['nrxns'])
        self.NSC_ci = np.zeros(self.runList[0].data.Reactions['nrxns'])
        for i in range(0, self.runList[0].data.Reactions['nrxns']):
            W = Wdata[:,2*i] + Wdata[:,2*i+1]             
            ci_info = Stats.cov_ci(W, TOFdata / self.TOF)
            self.NSC[i] = ci_info[0] + tof_fracs[2*i] + tof_fracs[2*i+1]
            self.NSC_ci[i] = ci_info[1]
                                   
    def WvarCheck(self): 
        
        ''' Compute trajectory derivative variances vs. time '''        
        
        W_dims = self.runList[0].data.Binary['W_sen_anal'].shape
        n_timepoints = W_dims[0]
        n_rxns = W_dims[1]       
        
        Wvars = np.zeros((n_timepoints,n_rxns))
        for i in range(0,n_timepoints):
            for j in range(0,n_rxns):
                data_vec = np.zeros((self.n_runs))
                for k in range(0,self.n_runs):
                    data_vec[k] = self.runList[k].data.Binary['W_sen_anal'][i,j]
                Wvars[i,j] = np.var(data_vec)
        
        ''' Plot results '''
        
        self.runList[0].PlotOptions()
        plt.figure()
            
        labels = []
        for i in range (2*len(self.runList[0].data.Reactions['names'])):
            if np.max(np.abs( Wvars[:,i] )) > 0:
                plt.plot(self.runList[0].data.Specnum['t'], Wvars[:,i])
                labels.append(self.runList[0].data.Reactions['names'][i/2])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('var(W)',size=24)
        plt.legend(labels,loc=4,prop={'size':20},frameon=False)        
        plt.show()
        
    def PlotSensitivities(self): 
        
        self.runAvg.PlotOptions()
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        
        for i in range (self.runList[0].data.Reactions['nrxns']):
            cutoff = 0.05
            if self.NSC[i] + self.NSC_ci[i] > cutoff or self.NSC[i] - self.NSC_ci[i] < -cutoff:     
                plt.barh(ind-0.9, self.NSC[i], width, color='r', xerr = self.NSC_ci[i], ecolor='k')               
                ylabels.append(self.runList[0].data.Reactions['names'][i])              
                yvals.append(ind-0.6)                
                ind = ind - 1

        plt.plot([0, 0], [0, ind], color='k')
        plt.xlim([0,1])
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('NSC',size=24)
        plt.yticks(yvals, ylabels)
        plt.show()
    
    def WriteSA_output(self,BatchPath):
        with open(BatchPath + 'SA_output.txt', 'w') as txt:
            txt.write('Normalized sensitivity coefficients \n\n')
            txt.write('Turnover frequency: \t' + '{0:.3E} \t'.format(self.TOF) + '+- {0:.3E} \t'.format(self.TOF_error) + '\n\n')               
            txt.write('Reaction name \t NSC \t NSC confidence \n')

            for rxn_ind in range(self.runList[0].data.Reactions['nrxns']):
                txt.write(self.runAvg.data.Reactions['names'][rxn_ind] + '\t' + '{0:.3f} +- \t'.format(self.NSC[rxn_ind]) + '{0:.3f}'.format(self.NSC_ci[rxn_ind]) + '\n')
                
    def FD_SA(self, rxn_inds = [1], pert_frac = 0.05, n_runs = 20, setup = True, exec_run = True, analyze_bool = True):
        
        # Create objects for perturbed systems
        plus = copy.deepcopy(self)
        minus = copy.deepcopy(self)
        FD_list = [plus, minus]

        for FD in FD_list:
            FD.n_runs = n_runs
        
        # Adjust pre-exponential factors in each
        adjust_plus = np.ones(self.runtemplate.data.Reactions['nrxns'])
        adjust_minus = np.ones(self.runtemplate.data.Reactions['nrxns'])
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