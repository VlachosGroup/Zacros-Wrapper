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
#Path is given as an example
Path = ut().SystemInformation()['Path']['Data'] + 'AtoB/'
Cnd = KMCut().InitializeCnd()
Cnd = RI().ReadAll(Path,Cnd)
Cnd = RO().ReadAll(Path,Cnd)
KMCut().pickleCnd(Cnd,'AtoB')




"""
This section of the code retreives and displays the python structure associated
with the Demo.p file
"""
#Cnd = KMCut().unpickleCnd('Demo')
#ut().PrintDict(Cnd)




"""
This section of the code perfoms a stiffness reduction on the Demo system for 
the given reduction modes and saves the result in a .p file
"""
TypeName = 'AtoB'
ModeList = ['linear_1_2','linear_2_4','tanh_1_2','tanh_2_4']

WallTimeList = [90]*len(ModeList)
for i in range(len(ModeList)):
    Mode = ModeList[i]
    WallTime = WallTimeList[i]

    if not os.path.isfile(ut().SystemInformation()['Path']['Data'] + 
        'PickledRunStructures/' + 'Recondition_' + TypeName +'_' + Mode + '.p'):
        RS().ReconditionCnd(CndIn = KMCut().unpickleCnd(TypeName),Name=TypeName+'_' + Mode 
                            ,RunParam = {'Event':1e3,'WallTime':WallTime,'Mode':Mode})




"""
This section builds jobs on Farber based on the prior stiffness reduction and
submits them
"""
nRareSample = 500 # Estimate, can be off depending on the system
nRuns = 2

for i in range(len(ModeList)):
    Mode = ModeList[i]
    OutList = KMCut().unpickleCnd('Recondition_' + TypeName+'_' + Mode)
    SDDict = KMCut().InitializeScaleDown()
    SDDict['SF'] = OutList['SF']
    SDDict['Mode'] = OutList['Cnd']['StiffnessRecondition']['Mode']
    SDDict['SDF'] = RS().CalculateScaleDown(BaseCnd='',Mode=SDDict['Mode'],SFIn=SDDict['SF'])['SDF']
    Cnd = OutList['Cnd']
    
    ETsS = np.log10(np.sum(SDDict['SF']/SDDict['SDF'])) # Estimated Timescale separation (for run time estimates)
    Cnd = KMCut().SetMaxEventNumber(Cnd,ETsS,nRareSample)
    JobPath = ut().SystemInformation()['Path']['Data'] + 'JobBuilds/' + TypeName+'_' + Mode + '/'
    if not os.path.isdir(JobPath):
        BI().BuildJob(Cnd,SDDict=SDDict,nRuns=nRuns,Name = TypeName+'_' + Mode)    
        if ut().SystemInformation()['OS'] == 'Linux':        
            p = subprocess.Popen("cd " + JobPath + "; " + "chmod 744 ./SubmitKMC.sh;./SubmitKMC.sh"
                            , stdout=subprocess.PIPE,shell=True)
            sys.stdout.flush()
        else:
            for j in ut().GetDir(JobPath):
                RunZacros().Run(JobPath + j + '/')




"""
This line goes through the job build folder, detects which jobs have completed
but have not been parsed, parses them, and saves the result
"""
RO().ReadAllJobs()


"""
This line reads the output of a parsed job which can then be processed.
"""
#CndList = RO().ReadJobOutput('Demo_linear_1_2')


"""
This section takes the various stiffness reduction runs and plots them to show
the effect of stiffness reduction on observed rate.
"""
nSites = 50
RxnStoich = [-1,1]
PropensityStoich = [-1,1,0,0]
LineColors = ['k','b']
#yLim = [0,0.0015]
RegFun = ['linear_','tanh_']
Modes = ['1_2','2_4']
PO().PlotReductionComparison(TypeName + '_',Modes,RegFun,nSites,RxnStoich,PropensityStoich,LineColors)#,yLim)
                            
                            
                            
                            
                            
                            
                            
                            
