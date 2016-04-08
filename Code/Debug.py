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
import numpy as np
import subprocess


#Path = 'E:/KMC_ToyProblem/'
#Path = ut().SystemInformation()['Path']['LocalRunDir'] + 'Run/'
#Cnd = KMCut().InitializeCnd()
#Cnd = RI().ReadAll(Path,Cnd)
#Cnd = RO().ReadAll(Path,Cnd)
#KMCut().pickleCnd(Cnd,'ToyProblem')

#Cnd = KMCut().unpickleCnd('Demo')
#ut().PrintDict(Cnd)

ModeList = ['tanh_1_2','tanh_2_4','tanh_3_6','tanh_4_8']
WallTimeList = [180,180,300,600]

for i in range(len(ModeList)):
    Mode = ModeList[i]
    WallTime = WallTimeList[i]

    RS().ReconditionCnd(CndIn = KMCut().unpickleCnd('ToyProblem'),Name='ToyProblem_' + Mode
                        ,RunParam = {'Event':5e3,'WallTime':WallTime,'Mode':Mode})


    OutList = KMCut().unpickleCnd('Recondition_ToyProblem_' + Mode)
    SDDict = KMCut().InitializeScaleDown()
    SDDict['SF'] = OutList['SF']
    SDDict['Mode'] = OutList['Cnd']['StiffnessRecondition']['Mode']
    SDDict['SDF'] = RS().CalculateScaleDown(BaseCnd='',Mode=SDDict['Mode'],SFIn=SDDict['SF'])['SDF']
    Cnd = OutList['Cnd']
    Cnd = KMCut().AdjustRuntime(Cnd,WallTime=3600*16)
    BI().BuildJob(Cnd,SDDict=SDDict,nRuns=10,Name = 'ToyProblem_' + Mode)
    if ut().SystemInformation()['OS'] == 'Linux':
        JobPath = ut().SystemInformation()['Path']['Data'] + 'JobBuilds/' + 'ToyProblem_' + Mode + '/'
        p = subprocess.Popen("cd " + JobPath + "; " + "chmod 744 ./SubmitKMC.sh;./SubmitKMC.sh"
                        , stdout=subprocess.PIPE,shell=True)



#JobPath = 'E:/ZacrosOutput/Runs/60Edge_MoreLateral_1atm/'
#CndList = RO().ReadJobOutput(JobPath)
#PO().CalcRate(CndList,[0.5,1.5,-1])
