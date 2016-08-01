# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 14:54:26 2016

@author: robieta
"""

from OutputData import OutputData
from RateRescaling import RateRescaling
import numpy as np

#import GeneralUtilities as ut
#import os
#import pickle
#import random
#import shutil
#import sys
#
#sys.path.append("..")
#from MachineSpecifics import MachineSpecifics as MS


class KMCrun:
    
    def __init__(self):
        
        self.output = OutputData()
        self.rescaling = RateRescaling()

    def PlotSpecVsTime(self):
        print 'Plot species vs. time'
        
    def PlotElemStepFreqs(self):
        print 'Plot elementary step frequencies'
    
    def ComputeTOF(self,Product):                       # return TOF and TOF error
        
        # Find the index of the product species
        product_ind = -1       
        for i in enumerate(self.output.input.Species['gas_spec']):
            if i[1] == Product:
                product_ind = i[0]
        
        # Make sure the index has been found
        if product_ind == -1:
            print 'Product species not found'
        else:
            product_ind = product_ind + self.output.input.Species['n_surf']         # Adjust index to account for surface species   
        
        
        nRxns = len(self.output.input.Reactions['Nu'])        
        TOF_contributions = [0 for i in range(nRxns)]              # number of product molecules produced in each reaction        
        for i, elem_stoich in enumerate(self.output.input.Reactions['Nu']):
            TOF_stoich = elem_stoich[product_ind]
            r = self.output.Binary['propCounter'][-1,i] / self.output.Specnum['t'][-1]
            TOF_contributions[i] = TOF_stoich * r

#            print TOF_stoich
#            print TOF            
            
        print TOF_contributions    
        TOF = np.sum(TOF_contributions)
        return TOF        
        
"""


      
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
            
"""