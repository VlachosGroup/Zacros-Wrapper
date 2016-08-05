# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 15:20:36 2016

@author: robieta
"""

from KMCrun import KMCrun
import GeneralUtilities as ut
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from Stats import Stats

##import ReadOutputFiles as RO
#import matplotlib.pyplot as plt

class AnalyzeData:
    
    def __init__(self):
                
        self.ParentFolder                     = ''
        self.runList                          = []      
        self.runAvg                           = KMCrun()            # Values are averages of all runs        
        
        self.ACF                              = {}
        self.ACF['Spacing']                   = {}
        self.ACF['Spacing']['Value']          = ''
        self.ACF['Spacing']['String']         = ''
        self.ACF['TauSep']                    = {}
        self.ACF['TauSep']['BurnIn']          = ''
        self.ACF['TauSep']['PostBurnIn']      = ''
        
        self.Product                          = ''
        self.TOF                              = ''
        self.TOF_error                        = ''
        self.NSC                              = ''
        self.NSC_ci                           = ''
        
        
    def ReadMultipleRuns(self,Path):
        print 'Reading output from ' + Path
        summary_fname = 'BatchSummary.p'
        if os.path.isfile(Path + summary_fname):                              # Runs have already been read
            print 'Existing pickle file found.'
            self.runList = pickle.load(open( Path + summary_fname, "rb" ))
        else:                                                               # Go through each directory and read the data           
            DirList = ut.GeneralUtilities().GetDir(Path)
            nDir = len(DirList)
            self.runList = ['' for i in range(nDir)]
            for i in range(nDir):
                print 'Reading run # ' + str(i+1) + ' / ' + str(nDir)
                RunPath = Path + DirList[i] + '/'
                self.runList[i]  = KMCrun()
                self.runList[i].output.Path = RunPath
                self.runList[i].output.ReadAllOutput()                               # Read input and output files
            pickle.dump( self.runList, open( Path + summary_fname, "wb" ) )          

    # Create a KMC run object with averaged species numbers, reaction firings, and propensities
    def AverageRuns(self):
        
        # Initialize run average with information from first run, then set data to zero
        self.runAvg = KMCrun()      
        self.runAvg.output.input = self.runList[0].output.input
        self.runAvg.output.Specnum['t'] = self.runList[0].output.Specnum['t']
             
        self.runAvg.output.Specnum['spec'] = self.runList[0].output.Specnum['spec'] - self.runList[0].output.Specnum['spec']         
        self.runAvg.output.Procstat['events'] = self.runList[0].output.Procstat['events'] - self.runList[0].output.Procstat['events']
        self.runAvg.output.Binary['cluster'] = self.runList[0].output.Binary['cluster'] - self.runList[0].output.Binary['cluster']
        self.runAvg.output.Binary['prop'] = self.runList[0].output.Binary['prop'] - self.runList[0].output.Binary['prop']
        self.runAvg.output.Binary['propCounter'] = self.runList[0].output.Binary['propCounter'] - self.runList[0].output.Binary['propCounter']        
        
        self.runAvg.output.Specnum['spec'] = self.runAvg.output.Specnum['spec'].astype(float)     
        self.runAvg.output.Procstat['events'] = self.runAvg.output.Procstat['events'].astype(float) 
        self.runAvg.output.Binary['cluster'] = self.runAvg.output.Binary['cluster'].astype(float)         
        
        # Add data from each run
        n_runs = len(self.runList)
        for run in self.runList:
            self.runAvg.output.Specnum['spec'] = self.runAvg.output.Specnum['spec'] + run.output.Specnum['spec'].astype(float) / n_runs         
            self.runAvg.output.Procstat['events'] = self.runAvg.output.Procstat['events'] + run.output.Procstat['events'].astype(float) / n_runs 
            self.runAvg.output.Binary['cluster'] = self.runAvg.output.Binary['cluster'] + run.output.Binary['cluster'].astype(float) / n_runs 
            self.runAvg.output.Binary['prop'] = self.runAvg.output.Binary['prop'] + run.output.Binary['prop'] / n_runs 
            self.runAvg.output.Binary['propCounter'] = self.runAvg.output.Binary['propCounter'] + run.output.Binary['propCounter'] / n_runs 

     
    def ComputeStats(self,product):
        self.TOF = self.runAvg.ComputeTOF(product)     
        
        n_runs = len(self.runList)
        n_rxns = len(self.runList[0].output.input.Reactions['Names'])
        Wdata = np.zeros((n_runs,n_rxns))      # number of runs x number of reactions
        TOFdata = np.zeros((n_runs))
        ind = 0
        for run in self.runList:
#            Wdata[ind,:] = run.output.Binary['W_sen_anal'][-1,:]
#            Wdata[ind,:] = run.output.Procstat['events'][-1,:] - run.output.Binary['propCounter'][-1,:]
            Wdata[ind,:] = run.output.Binary['W_sen_anal'][-1,:]            
            
                               
            
            TOFdata[ind] = run.ComputeTOF(product)
            ind = ind + 1             
        
        self.NSC = np.zeros((n_rxns/2,1))
        self.NSC_ci = np.zeros((n_rxns/2,1))
        for i in range(0,n_rxns/2):
            W = Wdata[:,2*i] + Wdata[:,2*i+1]
#            cov_mat = np.cov(W, TOFdata / self.TOF)     # normalize by the rate
#            self.NSC[i] = cov_mat[0,1]    
#            self.NSC_ci[i] = 0.1               
            ci_info = Stats().cov_ci(W, TOFdata / self.TOF)
            self.NSC[i] = ci_info[0]
            self.NSC_ci[i] = ci_info[1]
                                               
    def PlotSensitivities(self): 
        
        self.runAvg.PlotOptions()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        nrxns = len(self.runAvg.output.input.Reactions['Names'])
        for i in range (nrxns/2):
            cutoff = 0.05
            if self.NSC[i] + self.NSC_ci[i] > cutoff or self.NSC[i] - self.NSC_ci[i] < -cutoff:     
                plt.barh(ind-0.9, self.NSC[i], width, color='r', xerr = self.NSC_ci[i], ecolor='k')
                ylabels.append(self.runAvg.output.input.Reactions['Input'][i]['Name'])              
                yvals.append(ind-0.6)                
                ind = ind - 1

        plt.plot([0, 0], [0, ind], color='k')
        plt.xlim([-1 ,1])
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('NSC',size=24)
        plt.yticks(yvals, ylabels)
        plt.show()