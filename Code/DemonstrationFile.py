# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 14:00:23 2016

@author: robieta
"""

from BuildInputFiles import BuildInputFiles as BI
from KMCUtilities import KMCUtilities as KMCut
from ProcessOutput import ProcessOutput as PO
from ReadInputFiles import ReadInputFiles as RI
from ReadOutputFiles import ReadOutputFiles as RO
from ReconditionStiffness import ReduceStiffness as RS
from RunZacros import RunZacros
from GeneralUtilities import GeneralUtilities as ut
import sys
import os
import numpy as np
import subprocess
import pickle
import matplotlib.pyplot as plt


SystemName = 'WGS'
Mode = 'linear_1_2'
JobPath = ut().SystemInformation()['Path']['Data'] + 'JobBuilds/' + SystemName + '/'


"""
Read the files from a given path and writes the info to the Demo.p file for later use
"""
#Path = ut().SystemInformation()['Path']['Data'] + SystemName + '/'
#Cnd = KMCut().InitializeCnd()
#Cnd = RI().ReadAll(Path,Cnd)
#Cnd = RO().ReadAll(Path,Cnd)
#KMCut().pickleCnd(Cnd, SystemName)
#
#Cnd = KMCut().unpickleCnd(SystemName)       # retreive and displays the python structure associated with the Demo.p file     
#ut().PrintDict(Cnd)                         # print out information for the user to see

#print 'Starting reconditioning'        
#RS().ReconditionCnd(CndIn = KMCut().unpickleCnd(SystemName), Name = SystemName + '_unstiff',RunParam = {'Event':1e3,'MaxEvents': 100000,'Mode':Mode})


#"""
#Build jobs based on the prior stiffness reduction
#"""

#OutList = KMCut().unpickleCnd(SystemName + '_unstiff')
#Cnd = OutList['Cnd']
#
#SDDict = KMCut().InitializeScaleDown()
#SDDict['SF'] = OutList['SF']
#SDDict['Mode'] = OutList['Cnd']['StiffnessRecondition']['Mode']
#SDDict['SDF'] = OutList['SFList'][-1]
#
## Set conditions
#nRuns = 1000
#n_time_points = 100;                        # number of time sampled points to take
#Cnd['Conditions']['MaxStep'] = 'inf'
#Cnd['Conditions']['WallTime']['Max'] = 'inf'
#Cnd['Conditions']['SimTime']['Max'] = 0.5       # Need to choose this based off of the stiffness reduction
#Cnd['Report']['specnum'] = ['time', Cnd['Conditions']['SimTime']['Max'] / n_time_points]
#Cnd['StateInput']['Type'] = ''
#    
#BI().BuildJob(Cnd,SDDict=SDDict,nRuns=nRuns,Name = SystemName)
#
## Run all jobs (Windows only), alternatively transfer to Squidward and use job array
#for j in ut().GetDir(JobPath):
#    RunZacros().Run(JobPath + j + '/')

"""
Read KMC outputs
"""
#CndList = RO().ReadJobOutput(JobPath)                    # Use this if CndList does not exist yet

"""
Take the various stiffness reduction runs and plots them to show
the effect of stiffness reduction on observed rate.
"""

nSites = 1
#PropensityStoich = [0,-1,0,1]          # A-> B
PropensityStoich = [0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0]        # WGS

CndList = pickle.load(open( JobPath + 'CndList.p', "rb" ))
ut().PrintDict(CndList[0])
output = PO().CalcRateTransient(CndList, nSites, PropensityStoich)

print '\nReaction rate'
print output['Mean']

print '\nError'
print output['CI']

print '\nUnique Reaction stoichiometries'
print CndList[0]['Reactions']['UniqNu'][0::2]

print '\nNormalized sensitivity coefficients'
print output['SenCoeff']