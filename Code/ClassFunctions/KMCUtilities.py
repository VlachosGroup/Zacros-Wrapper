# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 14:54:26 2016

@author: robieta
"""

import copy
import GeneralUtilities as ut
import numpy as np
import os
import pickle
import random
from scipy import optimize
import shutil
import sys

sys.path.append("..")
from DynamicFiles.MachineSpecifics import MachineSpecifics as MS


class KMCUtilities:
    def __init__(self):
        pass
    
    def InitializeCnd(self,CndIn = {},PreserveInput = False, PreserveOutput = False):     
        Cnd = copy.copy(CndIn)
        if not PreserveInput:
            Cnd['Conditions'] = {}
            Cnd['Conditions']['T']                  = ''
            Cnd['Conditions']['P']                  = ''
            Cnd['Conditions']['Seed']               = ''
            Cnd['Conditions']['restart']            = ''
            Cnd['Conditions']['SimTime']            = {}
            Cnd['Conditions']['SimTime']['Max']     = ''
            Cnd['Conditions']['SimTime']['Actual']  = ''
            Cnd['Conditions']['WallTime']           = {}
            Cnd['Conditions']['WallTime']['Max']    = ''
            Cnd['Conditions']['WallTime']['Actual'] = ''
            Cnd['Conditions']['CPUTime']            = ''
            Cnd['Conditions']['nEvents']            = ''
            Cnd['Conditions']['MaxStep']            = ''
        
            Cnd['Species']                          = {}
            Cnd['Species']['n_gas']                 = ''
            Cnd['Species']['gas_spec']              = ''
            Cnd['Species']['gas_eng']               = ''
            Cnd['Species']['gas_MW']                = ''
            Cnd['Species']['gas_molfrac']           = ''
            Cnd['Species']['n_surf']                = ''
            Cnd['Species']['surf_spec']             = ''
            Cnd['Species']['surf_dent']             = ''
        
            Cnd['Report']                           = {}
            Cnd['Report']['specnum']                = ['','']
            Cnd['Report']['procstat']               = ['','']
            Cnd['Report']['hist']                   = ['','']
            Cnd['Report']['event']                  = ''
        
            Cnd['Cluster']                          = {}
            Cnd['Cluster']['nCluster']              = ''
            Cnd['Cluster']['nClusterVariant']       = ''
            Cnd['Cluster']['Input']                 = ''
            
            Cnd['Reactions']                        = {}
            Cnd['Reactions']['Names']               = ''
            Cnd['Reactions']['Nu']                  = ''
            Cnd['Reactions']['UniqNu']              = ''
            Cnd['Reactions']['Input']               = ''
            
            Cnd['StateInput']                       = {}
            Cnd['StateInput']['Type']               = ''
            Cnd['StateInput']['Struct']             = ''
            
            Cnd['Lattice']                          = {}
            Cnd['Lattice']['Input']                 = ''
            
            Cnd['StiffnessRecondition']             = {}
            Cnd['StiffnessRecondition']['Mode']     = ''
            Cnd['StiffnessRecondition']['APSdF']    = ''
        
        if not PreserveOutput:
            Cnd['Specnum']                          = {}
            Cnd['Specnum']['Spacing']               = ''
            Cnd['Specnum']['t']                     = ''
            Cnd['Specnum']['nEvents']               = ''
            Cnd['Specnum']['T']                     = ''
            Cnd['Specnum']['E']                     = ''
            Cnd['Specnum']['spec']                  = ''
            
            Cnd['Procstat']                         = {}
            Cnd['Procstat']['Spacing']              = ''
            Cnd['Procstat']['t']                    = ''
            Cnd['Procstat']['events']               = ''
            
            Cnd['History']                          = {}
            Cnd['History']['Final']                 = ''
            
            Cnd['Binary']                           = {}
            Cnd['Binary']['cluster']                = ''
            Cnd['Binary']['prop']                   = ''
            Cnd['Binary']['propCounter']            = ''
            
            Cnd['ACF']                              = {}
            Cnd['ACF']['Spacing']                   = {}
            Cnd['ACF']['Spacing']['Value']          = ''
            Cnd['ACF']['Spacing']['String']         = ''
            Cnd['ACF']['TauSep']                    = {}
            Cnd['ACF']['TauSep']['BurnIn']          = ''
            Cnd['ACF']['TauSep']['PostBurnIn']      = ''
  
        return Cnd
        
    def AdjustRuntime(self,Cnd,WallTime=24*3600):
        MaxLen = MS().MaxOutputEntries()
        EstimatedSimTime =  (Cnd['Conditions']['SimTime']['Actual']*WallTime/
                                Cnd['Conditions']['WallTime']['Actual'])
        Space1 = EstimatedSimTime / MaxLen
        Exponent = int(np.floor(np.log10(Space1)))
        Space2 = np.round(Space1 / 10 ** Exponent,0) * 10 ** Exponent
        
        Cnd['Conditions']['WallTime']['Max'] = WallTime
        Cnd['Conditions']['SimTime']['Max'] = ''
        
        Cnd['Report']['specnum']                = ['time',Space2]
        Cnd['Report']['procstat']               = ['time',Space2]
        Cnd['Report']['hist']                   = ['time',Space2*100]
        
        return Cnd
        
        
    def InitializeScaleDown(self):
        SDDict = {'SDF':'','SF':'','Mode':''}  
        return SDDict
        
    def KMCSeed(self):
        output = random.randint(1,int(2 ** 31 - 1))
        return output
        
    def pickleCnd(self,Cnd,Name):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        Path = sysinfo['Path']['Data'] + 'PickledRunStructures/'
        pickle.dump( Cnd, open( Path + Name + '.p', "wb" ) )
        
    def unpickleCnd(self,Name):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        Path = sysinfo['Path']['Data'] + 'PickledRunStructures/'
        if os.path.isfile(Path + Name + '.p'):
            Cnd = pickle.load(open( Path + Name + '.p', "rb" ))
        else:
            raise NameError('Specified pickle file does not exist')
        return Cnd
        
    def IsRun(self,Cnd,RequireBinaries = True):
         if Cnd['Specnum']['spec'] == '' or Cnd['Procstat']['events'] == '':
             RunBool = False
         elif RequireBinaries and Cnd['Binary'] != {}:
             BinBool1 = ut.GeneralUtilities().isblank(Cnd['Binary']['cluster'])
             BinBool2 = ut.GeneralUtilities().isblank(Cnd['Binary']['prop'])
             BinBool3 = ut.GeneralUtilities().isblank(Cnd['Binary']['propCounter'])
             if (BinBool1 or BinBool2 or BinBool3):
                 RunBool = False
             else:
                 RunBool = True
         else:
             RunBool = True
         return RunBool

    def CleanIntermediateRuns(self):
        RunDir = ut.GeneralUtilities().SystemInformation()['Path']['LocalRunDir']
        Folder = RunDir + 'IntermediateRuns/'
        if os.path.isdir(Folder):
            shutil.rmtree(Folder)
        os.mkdir(Folder)
        with open(Folder + 'Description.txt','w') as Txt:
            Txt.write('This folder is for storing intermediate results.\n')
            Txt.write('This folder is frequently purged. Do not rely on it for storing important data.')
        pass
    
    def CacheIntermediate(self,Name):
        Name = str(Name)
        RunDir = ut.GeneralUtilities().SystemInformation()['Path']['LocalRunDir']
        Files = ut.GeneralUtilities().GetFiles(RunDir + 'Run/')
        NewDir = RunDir + 'IntermediateRuns/' + Name + '/'
        if os.path.isdir(NewDir):
            shutil.rmtree(NewDir)
        os.mkdir(NewDir)
        for i in Files:
            shutil.copy(RunDir + 'Run/' + i,NewDir + i)
         
    def FindTau(self,InVec,tSpace=1,CoarseGrain = False):
        if np.min(InVec) == np.max(InVec): #Catch invariant species
            Tau = 0.
        else:
            Space = 1
            if CoarseGrain:
                nMin = MS().MinACFCoarseGrain()
                if len(InVec) > nMin:
                    Space = len(InVec)/nMin     # Use integer aritmetic
            AutoCorr = self.acf(InVec,Space=Space)
            xVec = np.array([float(i) * Space for i in range(0,len(AutoCorr))])
            SSError = lambda p,x,AutoCorr: sum((np.exp(-x/p)-AutoCorr
            ) ** 2)
            p1, success = optimize.leastsq(SSError, 1., args=(xVec,AutoCorr))
            Tau = float(p1 * tSpace)  
#            if Tau > 500:
#                pass
#                plt.plot(xVec,AutoCorr)

        return Tau
       
    def acf(self,InVec,lags=-1,Space=1):
        #InVec = InVec[None,:]
        if lags ==  -1 or lags > len(InVec)-2:
            lags = len(InVec)-2
#        AutoCorr = np.array([1.]*(lags))
        AutoCorr = []
        
        for i in range(1*Space, lags,Space):
            AutoCorrTemp = np.cov(np.append(InVec[None,:-i],InVec[None,i:],axis=0))[0,1]/np.var(InVec)
            AutoCorr.append(AutoCorrTemp)
        AutoCorr = np.array(AutoCorr)
        return AutoCorr
        
    def CndToTau(self,Cnd):
        if Cnd['Report']['specnum'][0] == 'time':
            if Cnd['Specnum']['spec'] != '':
                nCalc = 4
                SpecnumArray = np.array(Cnd['Specnum']['spec'])[:,:Cnd['Species']['n_surf']]
                TauList = [[[],[],[]] for i in range(nCalc)]
                MaxTauList = [['','',''] for i in range(nCalc)]
                StartInd = 0
                BurnInTauSepList = []
                for j in range(nCalc):
                    if j > 0:
                        StartInd = np.sum(BurnInTauSepList)
                        if SpecnumArray.shape[0] < StartInd:
                            break
                        
                    for i in range(SpecnumArray.shape[1]):
                        TauList[j][0].append(self.FindTau(SpecnumArray[StartInd:,i],1,CoarseGrain = True))
                    MaxTauList[j][0] = np.max(TauList[j][0])
                    
                    if Cnd['Binary']['cluster'] != '':
                        if j==0:
                            ClusterArray = np.array(Cnd['Binary']['cluster'])
                        for i in range(ClusterArray.shape[1]):
                            TauList[j][1].append(self.FindTau(ClusterArray[StartInd:,i],1,CoarseGrain = True))
                        MaxTauList[j][1] = np.max(TauList[j][1])
                     
                    if Cnd['Binary']['prop'] != '':
                        if j==0:
                            PropArray = np.array(Cnd['Binary']['prop'])
                        for i in range(PropArray.shape[1]):
                            TauList[j][2].append(self.FindTau(PropArray[StartInd:,i],1,CoarseGrain = True))
                        MaxTauList[j][2] = np.max(TauList[j][2])
                                       
                    if j < (nCalc-1):
                        BurnInTauSepList.append(int(np.ceil(np.max([i for i in MaxTauList[j] if i != '']))))
                        BurnInTau = np.max([i for i in MaxTauList[j] if i != ''])
    #                    ThreeTauBurnIn = int(np.ceil(BurnInTau*3))
                        ThreeTauBurnIn = int(np.ceil(BurnInTau))
                
                Cnd['ACF']['Spacing']['Value'] = Cnd['Specnum']['Spacing']
                if Cnd['ACF']['Spacing']['Value'] == 1:
                    StrNum = ''
                elif Cnd['ACF']['Spacing']['Value'] == 2:
                    StrNum = '2nd '
                elif Cnd['ACF']['Spacing']['Value'] == 3:
                    StrNum = '3rd '
                else:
                    StrNum = str(Cnd['ACF']['Spacing']['Value']) + 'th '
                Cnd['ACF']['Spacing']['String'] = 'Every ' + StrNum + 'value sampled'
                
                Cnd['ACF']['TauSep']['BurnIn'] =       [[-1 if i=='' else i for i in A] for A in MaxTauList[:-1]]
                Cnd['ACF']['TauSep']['PostBurnIn'] =   [-1 if i=='' else i for i in MaxTauList[-1]]

        return Cnd        
        
                
                
        