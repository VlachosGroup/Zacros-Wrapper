# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 14:00:23 2016

@author: robieta
"""

from ClassFunctions.BuildInputFiles import BuildInputFiles as BI
from ClassFunctions.KMCUtilities import KMCUtilities as KMCut
from ClassFunctions.ProcessOutput import ProcessOutput as PO
from ClassFunctions.ReadInputFiles import ReadInputFiles as RI
from ClassFunctions.ReadOutputFiles import ReadOutputFiles as RO
from ClassFunctions.ReconditionStiffness import ReduceStiffness as RS
from ClassFunctions.RunZacros import RunZacros
from ClassFunctions.GeneralUtilities import GeneralUtilities as ut
import sys
import os
import numpy as np
import subprocess
import pickle
import matplotlib.pyplot as plt

"""
This section of the code reads the files from a given path and writes the info
to the Demo.p file for later use
"""
#Path = ut().SystemInformation()['Path']['Data'] + 'NH3/'
#Cnd = KMCut().InitializeCnd()
#Cnd = RI().ReadAll(Path,Cnd)
#Cnd = RO().ReadAll(Path,Cnd)
#KMCut().pickleCnd(Cnd,'NH3')




"""
This section of the code retreives and displays the python structure associated
with the Demo.p file
"""
#Cnd = KMCut().unpickleCnd('AtoB')
#ut().PrintDict(Cnd)


"""
This section of the code perfoms a stiffness reduction on the Demo system for 
the given reduction modes and saves the result in a .p file
"""

#TypeName = 'AtoB'
#Mode = 'linear_1_2'
#MaxEvents = 10000
#
## See if the reconditioning has already been done
#if os.path.isfile(ut().SystemInformation()['Path']['Data'] + 'PickledRunStructures/' + 'Recondition_' + TypeName +'_' + Mode + '.p'): 
#    os.remove(ut().SystemInformation()['Path']['Data'] + 'PickledRunStructures/' + 'Recondition_' + TypeName +'_' + Mode + '.p')
#    
#print 'Starting reconditioning'        
#RS().ReconditionCnd(CndIn = KMCut().unpickleCnd(TypeName),Name=TypeName+'_' + Mode,RunParam = {'Event':1e3,'MaxEvents':MaxEvents,'Mode':Mode})



"""
This section builds jobs on Farber based on the prior stiffness reduction and
submits them
"""
#nRuns = 1000
#
#Mode = 'linear_1_2'
#OutList = KMCut().unpickleCnd('Recondition_' + TypeName+'_' + Mode)
#Cnd = OutList['Cnd']
#
#SDDict = KMCut().InitializeScaleDown()
#SDDict['SF'] = OutList['SF']
#SDDict['Mode'] = OutList['Cnd']['StiffnessRecondition']['Mode']
#SDDict['SDF'] = OutList['SFList'][-1]
#
#
## Set conditions
#Cnd['Conditions']['MaxStep'] = 'inf'
#Cnd['Conditions']['SimTime']['Max'] = 0.5
#Cnd['Conditions']['WallTime']['Max'] = 'inf'
#Cnd['Report']['specnum'] = ['time',0.005]
#Cnd['StateInput']['Type'] = ''
#    
#JobPath = ut().SystemInformation()['Path']['Data'] + 'JobBuilds/' + TypeName +'_' + Mode + '/'
#if not os.path.isdir(JobPath):
#    BI().BuildJob(Cnd,SDDict=SDDict,nRuns=nRuns,Name = TypeName+'_' + Mode)    
#    if ut().SystemInformation()['OS'] == 'Linux':        
#        p = subprocess.Popen("cd " + JobPath + "; " + "chmod 744 ./SubmitKMC.sh;./SubmitKMC.sh"
#                        , stdout=subprocess.PIPE,shell=True)
#        sys.stdout.flush()
#    else:
#        for j in ut().GetDir(JobPath):
#            RunZacros().Run(JobPath + j + '/')




"""
This line goes through the job build folder, detects which jobs have completed
but have not been parsed, parses them, and saves the result
"""
#RO().ReadAllJobs()


"""
This section takes the various stiffness reduction runs and plots them to show
the effect of stiffness reduction on observed rate.
"""
nSites = 1
PropensityStoich = [0,-1,0,1]
CndList = RO().ReadJobOutput('AtoB_linear_1_2')   
output = PO().CalcRateTransient(CndList, nSites, PropensityStoich)

print output['Mean']
print output['CI']
print output['SenCoeff']