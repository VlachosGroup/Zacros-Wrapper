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
from MachineSpecifics import MachineSpecifics as MS


class KMCUtilities:
    
    def __init__(self):
        
        self.Info = {}             # Dictionary class variable holds all simulation information       
        
        # --------------- Input data ---------------
        self.Info['Conditions']                       = {}
        self.Info['Conditions']['T']                  = ''
        self.Info['Conditions']['P']                  = ''
        self.Info['Conditions']['Seed']               = ''
        self.Info['Conditions']['restart']            = ''
        self.Info['Conditions']['SimTime']            = {}
        self.Info['Conditions']['SimTime']['Max']     = ''
        self.Info['Conditions']['SimTime']['Actual']  = ''
        self.Info['Conditions']['WallTime']           = {}
        self.Info['Conditions']['WallTime']['Max']    = ''
        self.Info['Conditions']['WallTime']['Actual'] = ''
        self.Info['Conditions']['CPUTime']            = ''
        self.Info['Conditions']['nEvents']            = ''
        self.Info['Conditions']['MaxStep']            = ''
    
        self.Info['Species']                          = {}
        self.Info['Species']['n_gas']                 = ''
        self.Info['Species']['gas_spec']              = ''
        self.Info['Species']['gas_eng']               = ''
        self.Info['Species']['gas_MW']                = ''
        self.Info['Species']['gas_molfrac']           = ''
        self.Info['Species']['n_surf']                = ''
        self.Info['Species']['surf_spec']             = ''
        self.Info['Species']['surf_dent']             = ''
    
        self.Info['Report']                           = {}
        self.Info['Report']['specnum']                = ['','']
        self.Info['Report']['procstat']               = ['','']
        self.Info['Report']['hist']                   = ['','']
        self.Info['Report']['event']                  = ''
    
        self.Info['Cluster']                          = {}
        self.Info['Cluster']['nCluster']              = ''
        self.Info['Cluster']['nClusterVariant']       = ''
        self.Info['Cluster']['Input']                 = ''
        
        self.Info['Reactions']                        = {}
        self.Info['Reactions']['Names']               = ''
        self.Info['Reactions']['Nu']                  = ''
        self.Info['Reactions']['UniqNu']              = ''
        self.Info['Reactions']['Input']               = ''
        
        self.Info['StateInput']                       = {}
        self.Info['StateInput']['Type']               = ''
        self.Info['StateInput']['Struct']             = ''
        
        self.Info['Lattice']                          = {}
        self.Info['Lattice']['Input']                 = ''
        
        self.Info['StiffnessRecondition']             = {}
        self.Info['StiffnessRecondition']['Mode']     = ''
        self.Info['StiffnessRecondition']['APSdF']    = ''
        
        # --------------- Output data ---------------
        self.Info['Specnum']                          = {}
        self.Info['Specnum']['Spacing']               = ''
        self.Info['Specnum']['t']                     = ''
        self.Info['Specnum']['nEvents']               = ''
        self.Info['Specnum']['T']                     = ''
        self.Info['Specnum']['E']                     = ''
        self.Info['Specnum']['spec']                  = ''
        
        self.Info['Procstat']                         = {}
        self.Info['Procstat']['Spacing']              = ''
        self.Info['Procstat']['t']                    = ''
        self.Info['Procstat']['events']               = ''
        
        self.Info['History']                          = {}
        self.Info['History']['Final']                 = ''
        
        self.Info['Binary']                           = {}
        self.Info['Binary']['cluster']                = ''
        self.Info['Binary']['prop']                   = ''
        self.Info['Binary']['propCounter']            = ''
        self.Info['Binary']['W_sen_anal']             = ''
        
        self.Info['ACF']                              = {}
        self.Info['ACF']['Spacing']                   = {}
        self.Info['ACF']['Spacing']['Value']          = ''
        self.Info['ACF']['Spacing']['String']         = ''
        self.Info['ACF']['TauSep']                    = {}
        self.Info['ACF']['TauSep']['BurnIn']          = ''
        self.Info['ACF']['TauSep']['PostBurnIn']      = ''
  
        # --------------- Analysis data ---------------  

        
    def BuildEmptyFolders(self):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        if not os.path.isdir(sysinfo['Path']['LocalRunDir']):
            os.mkdir(sysinfo['Path']['LocalRunDir'])
            os.mkdir(sysinfo['Path']['LocalRunDir'] + 'Build/')
            os.mkdir(sysinfo['Path']['LocalRunDir'] + 'Run/')
            os.mkdir(sysinfo['Path']['LocalRunDir'] + 'IntermediateRuns/')
        
        
    def AdjustRuntime(self,WallTime=24*3600):
        MaxLen = MS().MaxOutputEntries()
        EstimatedSimTime =  (self['Conditions']['SimTime']['Actual']*WallTime/
                                self['Conditions']['WallTime']['Actual'])
        Space1 = EstimatedSimTime / MaxLen
        Exponent = int(np.floor(np.log10(Space1)))
        Space2 = np.round(Space1 / 10 ** Exponent,0) * 10 ** Exponent
        
        self.Info['Conditions']['WallTime']['Max'] = WallTime
        self.Info['Conditions']['SimTime']['Max'] = ''
        
        self.Info['Report']['specnum']                = ['time',Space2]
        self.Info['Report']['procstat']               = ['time',Space2]
        self.Info['Report']['hist']                   = ['time',Space2*100]
        
        return self
        
    def SetMaxEventNumber(self,ETsS,nRare):
        # ETsS = Estimated Timescale Separation (log10)
        MaxLen = MS().MaxOutputEntries()
        MaxEvents = (10 ** ETsS) * nRare
        self['Conditions']['MaxStep'] = int(MaxEvents)
        self['Conditions']['SimTime']['Max'] = 'inf'
        self['Conditions']['WallTime']['Max'] = 'inf'
        EstimatedSimTime = (self['Conditions']['SimTime']['Actual']
                            /self['Conditions']['nEvents']*MaxEvents)
        Space1 = EstimatedSimTime / MaxLen
        Exponent = int(np.floor(np.log10(Space1)))
        Space2 = np.round(Space1 / 10 ** Exponent,0) * 10 ** Exponent
        
        self['Report']['specnum']                = ['time',Space2]
        self['Report']['procstat']               = ['time',Space2]
        self['Report']['hist']                   = ['time',Space2*100]
        
        return self
        
    
    def InitializeScaleDown(self):
        SDDict = {'SDF':'','SF':'','Mode':''}  
        return SDDict
        
    def KMCSeed(self):
        output = random.randint(1,int(2 ** 31 - 1))
        return output
        
    def pickleself(self,Name):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        Path = sysinfo['Path']['Data'] + 'Pickles/'
        pickle.dump( self, open( Path + Name + '.p', "wb" ) )
        
    def unpickleself(self,Name):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        Path = sysinfo['Path']['Data'] + 'Pickles/'
        print Path + Name + '.p'
        if os.path.isfile(Path + Name + '.p'):
            self = pickle.load(open( Path + Name + '.p', "rb" ))
        else:
            raise NameError('Specified pickle file does not exist')
        return self
        
    def IsRun(self,RequireBinaries = True):
         if self['Specnum']['spec'] == '' or self['Procstat']['events'] == '':
             RunBool = False
         elif RequireBinaries and self['Binary'] != {}:
             BinBool1 = ut.GeneralUtilities().isblank(self['Binary']['cluster'])
             BinBool2 = ut.GeneralUtilities().isblank(self['Binary']['prop'])
             BinBool3 = ut.GeneralUtilities().isblank(self['Binary']['propCounter'])
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
        self.BuildEmptyFolders()
        Name = str(Name)
        RunDir = ut.GeneralUtilities().SystemInformation()['Path']['LocalRunDir']
        Files = ut.GeneralUtilities().GetFiles(RunDir + 'Run/')
        NewDir = RunDir + 'IntermediateRuns/' + Name + '/'
        if os.path.isdir(NewDir):
            shutil.rmtree(NewDir)
        os.mkdir(NewDir)
        for i in Files:
            shutil.copy(RunDir + 'Run/' + i,NewDir + i)