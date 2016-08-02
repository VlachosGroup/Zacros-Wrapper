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
        self.NSC_error                        = ''
        
        
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
     
        n_rxns = len(self.runList[0].output.input.Reactions['Names'])
        for i in range(0,n_rxns):
            print ' '
            print self.runList[0].output.input.Reactions['Names'][i]
            print self.runAvg.output.Procstat['events'][-1,i]
            print self.runAvg.output.Binary['propCounter'][-1,i] 
     
    def ComputeStats(self,product):
        self.TOF = self.runAvg.ComputeTOF(product)
        
        print 'TOF'
        print self.TOF        
        
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
        
        print ' '
        print 'NSC'        
        
        NSC = np.zeros((n_rxns/2,1))
        NSC_ci = np.zeros((n_rxns/2,1))
        for i in range(0,n_rxns/2):
            W = Wdata[:,2*i] + Wdata[:,2*i+1]
            cov_mat = np.cov(W, TOFdata)
            NSC[i] = cov_mat[0,1] / self.TOF                   # normalize by the rate
            NSC_ci[i] = 0                                   # need to implement statistical bootstrapping to estimate the confidence interval
            
            print ' '            
            print self.runList[0].output.input.Reactions['Names'][2*i]
            print NSC[i]