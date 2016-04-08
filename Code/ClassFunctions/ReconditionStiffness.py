# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 20:28:48 2016

@author: RDX
"""

import BuildInputFiles as BI
import copy
import GeneralUtilities as ut
import KMCUtilities as KMCut
import numpy as np
import os
import ReadInputFiles as RI
import ReadOutputFiles as RO
import RunZacros as RunZacros

class ReduceStiffness:
    def __init__(self):
        pass
    
    def DefaultRunParam(self):
        RunParam = {'Event':1e4,'WallTime':60,'Mode':'tanh_2_4'}
        return RunParam
        
    def ReconditionCnd(self,CndIn,RunParam = '',nSim = 50,Name=''):
        RunBool = KMCut.KMCUtilities().IsRun(CndIn)

        RunPath = ut.GeneralUtilities().SystemInformation()['Path']['LocalRunDir'] + 'Run/'
        CndList = []
        if RunParam == '':
            RunParam = self.DefaultRunParam()
        
        
        ReportCheck1 = (CndIn['Report']['specnum'][0] == 'event' and 
                        CndIn['Report']['procstat'][0] == 'event' and 
                        CndIn['Report']['hist'][0] == 'event')
                        
        ReportCheck2 = (CndIn['Report']['specnum'][1] == CndIn['Report']['procstat'][1] and
                        CndIn['Report']['specnum'][1] == CndIn['Report']['hist'][1]
                        and CndIn['Report']['specnum'][1] == RunParam['WallTime'])
                        
        ReportCheck3 = not ut.GeneralUtilities().isblank(CndIn['History']['Final'])
                      
        if not (ReportCheck1 and ReportCheck2 and ReportCheck3):
            CndIn['Report']['specnum'][0] = 'event'
            CndIn['Report']['specnum'][1] = RunParam['Event']
            CndIn['Report']['procstat'][0] = 'event'
            CndIn['Report']['procstat'][1] = RunParam['Event']
            CndIn['Report']['hist'][0] = 'event'
            CndIn['Report']['hist'][1] = RunParam['Event']
            CndIn = KMCut.KMCUtilities().InitializeCnd(CndIn,PreserveInput = True)
            RunBool = False               
        
        if not RunBool:
            CndIn['Conditions']['WallTime']['Max'] = RunParam['WallTime']
            BI.BuildInputFiles().BuildFiles(CndIn)
            RunZacros.RunZacros().Run()
            CndIn = RO.ReadOutputFiles().ReadAll(RunPath,CndIn)
            
        CndList.append(CndIn)
        
        CndClean = KMCut.KMCUtilities().InitializeCnd(CndIn,PreserveInput = True)    
        SFList = []
        KMCut.KMCUtilities().CleanIntermediateRuns()
        Converged = False
        for i in range(0,nSim):
            SDDict = self.CalculateScaleDown(CndList[i],RunParam['Mode'])
            Cnd = copy.copy(CndClean)
            Cnd['StateInput']['Type'] = 'history'
            Cnd['StateInput']['Struct'] = CndList[i]['History']['Final']
            Cnd['Conditions']['Seed'] = ''            
            BI.BuildInputFiles().BuildFiles(Cnd,SDDict = SDDict)
            RunZacros.RunZacros().Run()
            KMCut.KMCUtilities().CacheIntermediate(i+1)
            Cnd = RI.ReadInputFiles().ReadAll(RunPath)
            Cnd = RO.ReadOutputFiles().ReadAll(RunPath,Cnd)
            CndList.append(copy.copy(Cnd))
            SFList.append(self.CalculateScaleDown(CndList[i+1],RunParam['Mode'])['SF'])
            print SFList[-1]
            Converged,SFOut = self.TestConvergence(SFList)
            print 'Pass ' + str(i)
            if Converged:
                PickleName = 'Recondition'
                if Name != '':
                    PickleName += '_' + str(Name)
                KMCut.KMCUtilities().pickleCnd({'Cnd':Cnd,'SF':SFOut},PickleName)
                break

    def TestConvergence(self,SFList):
        Converged = False
        MinPostSample = 5
        SFOut = ''

        nRuns = len(SFList)
        SF = np.array(copy.copy(SFList))
        nRxn = SF.shape[1]
        TauMax = 0.
        for i in range(nRxn):
            Tau = KMCut.KMCUtilities().FindTau(np.log10(SF[:,i]))
            if not np.isnan(Tau):
                if Tau > TauMax:
                    TauMax = copy.copy(Tau)
        StartInd = int(np.ceil(TauMax*6))

        if nRuns > StartInd:
            TauMax2 = 0.
            for i in range(nRxn):
                Tau = KMCut.KMCUtilities().FindTau(np.log10(SF[StartInd:,i]))
                if not np.isnan(Tau):
                    if Tau > TauMax2:
                        TauMax2 = copy.copy(Tau)
            MinPostSample = np.max([np.ceil(TauMax2*3),MinPostSample])
            Converged = (nRuns - StartInd) > MinPostSample
            if Converged:
                LSF = np.log10(np.array([[np.max([i,1]) for i in SF[j,:].tolist()] for j in range(StartInd,nRuns)]))
                MLSF = np.mean(LSF,axis=0).tolist()
                LSFCI = ut.GeneralUtilities().CI(LSF,axis=0).tolist()
                MaxCIOverMean = 0
                for i in range(len(MLSF)):
                    if MLSF[i] > 1:
                        MaxCIOverMean = np.max([MaxCIOverMean,LSFCI[i]/MLSF[i]])

                if MaxCIOverMean > 1:
                    Converged = False
                    
        if Converged:
            SFOut = (10 ** np.array(MLSF)).tolist()
        return Converged,SFOut

        

    
    def CalculateScaleDown(self,BaseCnd,Mode,SFIn = ''):
        TransformFunction = Mode.split('_')[0]
        Cutoff = float(Mode.split('_')[1])
        UpperLimit = float(Mode.split('_')[2])
        
        if SFIn == '':        
            Species = BaseCnd['Species']['surf_spec'] + BaseCnd['Species']['gas_spec']
            nEvents = np.sum(BaseCnd['Procstat']['events'][-1])
            
            UniqNu = np.array(BaseCnd['Reactions']['UniqNu'])
            for i in range(UniqNu.shape[0]-1,0,-1):
                if np.array_equal(UniqNu[i],-UniqNu[i-1]):
                    UniqNu = np.delete(UniqNu,(i),axis=0)
               
            MechDict = BaseCnd['Reactions']['Input']
            Nu = np.array([[0]*len(Species)]*len(MechDict))
            NuInd = []       
            Count = -1
    
            for i in MechDict:
                Count += 1
                Count2 = 1
                for j in i['initial']:
                    if int(j.split()[0]) == Count2:
                        Count2 += 1
                        for k in range(0,len(Species)):
                            if j.split()[1] == Species[k]:
                                Nu[Count,k] += -1
                Count2 = 1
                for j in i['final']:
                    if int(j.split()[0]) == Count2:
                        Count2 += 1
                        for k in range(0,len(Species)):
                            if j.split()[1] == Species[k]:
                                Nu[Count,k] += 1
                                
                if len(i['gas_reacs_prods']) > 0:
                    for k in range(0,len(Species)):
                        if i['gas_reacs_prods'][0] == Species[k]:
                            Nu[Count,k] += int(i['gas_reacs_prods'][1])          
                
                for j in range(0,UniqNu.shape[0]):
                    if np.array_equal(Nu[Count,:],UniqNu[j]) or np.array_equal(Nu[Count,:],-UniqNu[j]):
                        for k in range(0,len(i['variant'])):
                            NuInd.append(j)
            
            if BaseCnd['StiffnessRecondition']['APSdF'] == '':
                APSdF = np.array([1. for i in range(0,BaseCnd['Binary']['propCounter'].shape[1]/2)])
            else:
                APSdF = np.array(BaseCnd['StiffnessRecondition']['APSdF'])
                
            NuGroups = np.unique(NuInd)
            GroupProp_ScaleupCorrected = np.array([[0.]*2]*len(NuGroups))
            Count = -1
            for i in NuGroups:
                Count += 1
                Count2 = -1
                for j in NuInd:
                    Count2 += 1
                    if i==j:
                        GroupProp_ScaleupCorrected[Count,0] += (
                            BaseCnd['Binary']['propCounter'][-1,0::2][Count2] * APSdF[Count2])
                        GroupProp_ScaleupCorrected[Count,1] += (
                            BaseCnd['Binary']['propCounter'][-1,1::2][Count2] * APSdF[Count2])
                            
            ApproachToEquil = 1 - (np.abs(BaseCnd['Binary']['propCounter'][-1,0::2]
                                    - BaseCnd['Binary']['propCounter'][-1,1::2])/
                                    (BaseCnd['Binary']['propCounter'][-1,0::2] 
                                    + BaseCnd['Binary']['propCounter'][-1,1::2]+1./nEvents/APSdF))
    
            IntRxnSpeed = np.max(GroupProp_ScaleupCorrected,axis=1)
            PropSlowScale = np.max(IntRxnSpeed)     
            for i in IntRxnSpeed:
                if i > 0. and i < PropSlowScale:
                    PropSlowScale = i
                    
            DegOfSep = (np.max(np.array([BaseCnd['Binary']['propCounter'][-1,0::2],
                         BaseCnd['Binary']['propCounter'][-1,1::2]]),axis=0)
                        * APSdF / PropSlowScale)
                        
            StiffnessFactor = (ApproachToEquil ** 2.) * DegOfSep
        else:
            StiffnessFactor = np.array(SFIn)
            APSdF = [1. for i in range(StiffnessFactor.shape[0])]
            nEvents = np.max(StiffnessFactor) ** 2
        
        StiffInd = []
        for i in range(0,StiffnessFactor.shape[0]):
            if StiffnessFactor[i] >  10 ** Cutoff or APSdF[i] > 1.:
                StiffInd.append(i)

        SDF_out = np.array([1.]*len(StiffnessFactor))
        SF = np.array(StiffnessFactor[StiffInd])
        LSF = np.log10(SF)      # Log stiffness factor
        ELSF = LSF - Cutoff     # Excess log stiffness factor
        
        # Regularize ELSF
        if TransformFunction == 'linear':
            RELSF = ELSF / np.max(ELSF)
        elif TransformFunction == 'tanh':
            RELSF = np.tanh(np.e*ELSF/UpperLimit) / np.tanh(np.e*np.max(ELSF)/UpperLimit)
            
        if (UpperLimit - Cutoff) < np.max(ELSF):
            TELSF = (UpperLimit - Cutoff) * RELSF   # Transform excess log stiffness factor
            TLSF = TELSF + Cutoff                   # Transform log stiffness factor
            SDF = SF / (10**TLSF)                   # Scaledown factor
            for i in range(0,len(StiffInd)):
                if SDF[i]/APSdF[StiffInd[i]] > np.sqrt(nEvents) / 2.:  # Prevent Runaway scaledown
                    SDF[i] = APSdF[StiffInd[i]] * np.sqrt(nEvents) / 2.
                elif SDF[i]/APSdF[StiffInd[i]] < 1./np.sqrt(nEvents) * 2.:  # Prevent Runaway scaledown
                    SDF[i] = APSdF[StiffInd[i]] / np.sqrt(nEvents) * 2.
            SDF_out[StiffInd] = SDF
            
        for i in range(0,len(SDF_out)):
            if SDF_out[i] < 1.0:
                SDF_out[i] = 1.0    
            
        SDDict = {'SDF':SDF_out,'SF':StiffnessFactor,'Mode':Mode}
        return SDDict
        






