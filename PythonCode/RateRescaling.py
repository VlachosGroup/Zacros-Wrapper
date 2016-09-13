# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 20:28:48 2016

@author: RDX
"""

#import BuildInputFiles as BI
#import copy
#import GeneralUtilities as ut
#import numpy as np
#import os
#import ReadInputFiles as RI
#import ReadOutputFiles as RO
#import RunZacros as RunZacros

import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np
from KMCrun import KMCrun

class RateRescaling:
    
    def __init__(self):
        
        self.KMC_system = KMCrun()
        self.SDF_mat    = []        # scaledown factors for each iteration
    
    def PerformScaledown(self):
        # Print reaction names and scaledown factors into a file
        print 'Rescaling rate constants\n'
        
        max_iterations = 15
        max_events = 1e5
        converged = False
        iteration = 0        
        cutoff = 0.5                # If a rate constant is changing by a half-order of magnititude or more, continue        
        
        self.SDF_mat = self.KMC_system.data.scaledown_factors
        
        while not converged and iteration < max_iterations:
            print 'Iteration number ' + str(iteration)
            print self.KMC_system.data.scaledown_factors
            print '\n'
            iteration += 1
            
            # Run KMC simulation
            delta_sdf = self.ProcessStepFreqs()         # compute change in scaledown factors based on simulation result
            # Check convergence
#            if np.max(np.log10(delta_sdf)) < cutoff:             # if delta_sdf's are small
#                converged = True
            
            self.SDF_mat = np.vstack([self.SDF_mat,self.KMC_system.data.scaledown_factors])            
            for i in range (len(self.KMC_system.data.Reactions['nrxns'])):            
            
            
        self.PlotStiffnessReduction()
        print self.SDF_mat
        # return time-scale information, use this to set a good sampling time for big run
                  
    def ProcessStepFreqs(self):
        # Process KMC output and determine how to further scale down reactions
        print 'Rescaling rate constants\n'
        self.KMC_system.data.ReadAllOutput()
        delta_sdf = np.ones(self.KMC_system.data.Reactions['nrxns'])
        return delta_sdf
 
    def WriteScaledownSummary(self,flname):
        # Print reaction names and scaledown factors into a file
        print 'Write scaledown summary'        
        
    def ReadScaledownSummary(self,flname):
        # Read reaction names and scaledown factors from a file
        print 'Read scaledown summary'
    
    def PlotStiffnessReduction(self):
        
        ''' Data '''
        iterations = np.array([0, 1, 2, 3, 4, 5])
        COads_SDF = np.array([1, 1.1e-2, 8.9e-4, 5.0e-4, 3.3e-4, 3.1e-4])
        H2Oads_SDF = np.array([1, 8.2e-3, 1.1e-4, 9.3e-6, 3.3e-6, 2.1e-6])
        slow_rxn = np.array([1, 1, 1, 1, 1, 1])
        
        ''' Plotting '''
        
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 16
        
        plt.figure()
        
        plt.plot(iterations, COads_SDF, 'o-', markersize = 15)
        plt.plot(iterations, H2Oads_SDF, 'o-', markersize = 15)
        plt.plot(iterations, slow_rxn, 'o-', markersize = 15)
        
        plt.xticks(size=24)
        plt.yticks([1e-6, 1e-4, 1e-2, 1e0], size=24)            # not working
        plt.xlabel('iterations',size=30)
        plt.ylabel('scaledown factor',size=30)
        plt.legend(['CO ads.', '$H_2O$  $ads.$', 'slow rxns'],loc=1,prop={'size':20},frameon=False)
        plt.show()
        
        plt.yscale('log')
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)