# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 15:20:36 2016

@author: robieta
"""

from KMCrun import KMCrun
import GeneralUtilities as ut
import os
import pickle

##import ReadOutputFiles as RO
#import matplotlib.pyplot as plt
#import numpy as np

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

    def AverageRuns(self):
        self.runAvg = self.runList[0]       # Placeholder: Really want to average all of the runs, not just copy the first one
        
"""   
        
    # Should perhaps make the sensitivity analysis a separate function so we can use it for both the transient and steady-state analysis    
    def CalcRateTransient(self,CndIn,nSites=1,PropStoich=''):
        if type(CndIn) == dict:
            CndIn = [CndIn]
            
        n_rxns = len(CndIn[0]['Reactions']['Names'])
        WList = np.zeros((len(CndIn),n_rxns))      # number of runs x number of reactions
        rate_contributions = np.zeros((len(CndIn),n_rxns))      # number of runs x number of reactions
        ind = 0

        for Cnd in CndIn:        
                
            if PropStoich != '':               
                
                for i,PS in enumerate(Cnd['Reactions']['Nu']):
                    if np.array_equal(PS,PropStoich):
                        rate_contributions[ind,i] = Cnd['Binary']['prop'][-1,i]
                    elif np.array_equal(PS,-np.array(PropStoich)):
                        rate_contributions[ind,i] = -Cnd['Binary']['prop'][-1,i]
                    else:
                        rate_contributions[ind,i] = 0
            else:
                print 'Reaction rate stoichiometry not specified'                        
            WList[ind,:] = np.array(Cnd['Binary']['W_sen_anal'][-1,:])      # Store sesitivity analysis data
            ind += 1
            
        avg_rates = np.mean(rate_contributions,axis=0)
        rate_fracs = avg_rates / np.sum(avg_rates)
        rate_obj = np.sum(rate_contributions,axis=1)        
        
        # Compute averages and confidence intervals
        rate_obj = rate_obj / nSites            # normalize the observable by the number of active sites
        Mean = np.mean(rate_obj)
        CI = ut.GeneralUtilities().CI(np.array(rate_obj),axis=0)          
        
        # Compute sensitivity coefficients
        SenCoeffraw = np.zeros(n_rxns)
        SenCoeffCI = np.zeros(n_rxns)
        
        for i in range(0,n_rxns):
#            if CndIn[0]['StiffnessRecondition']['APSdF'][int(np.floor(i/2))] > 1:       # Reaction is fast and unimportant 
#                SenCoeffraw[i] = 0
#                SenCoeffCI[i] = 0
#            else:                                           # Reaction is slow and may be important
#                W = WList[:,i]                              # need to group these by reaction class
#                cov_mat = np.cov(W,rate_obj)
#                SenCoeffraw[i] = cov_mat[0,1] / Mean                   # normalize by the rate
#                SenCoeffCI[i] = 0                                   # need to implement statistical bootstrapping to estimate the confidence interval
        
            W = WList[:,i]                              # need to group these by reaction class
            cov_mat = np.cov(W,rate_obj)
            SenCoeffraw[i] = cov_mat[0,1] / Mean                   # normalize by the rate
            SenCoeffCI[i] = 0                                   # need to implement statistical bootstrapping to estimate the confidence interval        
        
        SenCoeffraw = SenCoeffraw + rate_fracs    # Add extra term which accounts for changes in the rate constant changing the propensity value      
        
        # Combine NSCs for like reactions
        
        unique_rxns = CndIn[0]['Reactions']['UniqNu']
        SenCoeffCombined = np.zeros(int(len(unique_rxns)/2))
        for i in range(0,len(SenCoeffCombined)):
            for j,PS in enumerate(CndIn[0]['Reactions']['Nu']):
                if np.array_equal(PS,unique_rxns[2*i]) | np.array_equal(PS,unique_rxns[2*i+1]):
                    SenCoeffCombined[i] = SenCoeffCombined[i] + SenCoeffraw[j]
                        
#            print unique_rxns[2*i]
#            print unique_rxns[2*i+1]
#            print 'hi'
            # unique reactions are in a weird order...
            
        Output = {'Mean':Mean,'CI':CI,'SenCoeff':SenCoeffCombined,'SenCoeffCI':SenCoeffCI}

        return Output
"""