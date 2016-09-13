# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 14:54:26 2016

@author: robieta
"""

from KMCrun_data import KMCrun_data
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat

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
        
        self.data = KMCrun_data()

    def PlotOptions(self):
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 12

    def PlotSurfSpecVsTime(self):
        self.PlotOptions()
        plt.figure()
        
        for i in range (len(self.data.Species['surf_spec'])):
            plt.plot(self.data.Specnum['t'], self.data.Specnum['spec'][:,i])    
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('spec. pop.',size=24)
#        plt.ylabel('coverage',size=24)
        plt.legend(self.data.Species['surf_spec'],loc=4,prop={'size':20},frameon=False)        
        plt.show()
    
    def PlotGasSpecVsTime(self):
        self.PlotOptions()
        plt.figure()          
          
        for i in range (len(self.data.Species['gas_spec'])):
            ind = i + len(self.data.Species['surf_spec'])
            plt.plot(self.data.Specnum['t'], self.data.Specnum['spec'][:,ind])    
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('spec. pop.',size=24)
        plt.legend(self.data.Species['gas_spec'],loc=2,prop={'size':20},frameon=False)        
        plt.show()     
    
    def PlotPropsVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.data.Reactions['Names'])):
            if np.max(np.abs(self.data.Binary['prop'][:,i])) > 0:
                plt.plot(self.data.Specnum['t'], self.data.Binary['prop'][:,i]) 
                labels.append(self.data.Reactions['Names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('props',size=24)
        plt.legend(self.data.Reactions['Names'],loc=2,prop={'size':20},frameon=False)        
        plt.show()
        
    def PlotIntPropsVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.data.Reactions['Names'])):
            if np.max(np.abs(self.data.Binary['propCounter'][:,i])) > 0:
                plt.plot(self.data.Specnum['t'], self.data.Binary['propCounter'][:,i]) 
                labels.append(self.data.Reactions['Names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('integral props',size=24)
        plt.legend(self.data.Reactions['Names'],loc=2,prop={'size':20},frameon=False)        
        plt.show()
    
    def PlotWVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.data.Reactions['Names'])):
            if np.max(np.abs(self.data.Binary['W_sen_anal'][:,i])) > 0:
                plt.plot(self.data.Specnum['t'], self.data.Binary['W_sen_anal'][:,i]) 
#                plt.plot(self.data.Specnum['t'], self.data.Procstat['events'][:,i] - self.data.Binary['propCounter'][:,i]) 
                labels.append(self.data.Reactions['Names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('W',size=24)
        plt.legend(self.data.Reactions['Names'],loc=2,prop={'size':20},frameon=False)        
        plt.show()       
    
    def PlotElemStepFreqs(self):
        self.PlotOptions
        plt.figure()        
        
        width = 0.2
        ind = 0
        yvals = []
        ylabels = []
        nRnxs = len(self.data.Reactions['Names'])

        for i in range (nRnxs/2): 
            if self.data.Procstat['events'][-1,2*i] + self.data.Procstat['events'][-1,2*i+1] > 0:
                net_freq = abs(self.data.Procstat['events'][-1,2*i] - self.data.Procstat['events'][-1,2*i+1])               
                if self.data.Procstat['events'][-1,2*i] > 0:              
                    plt.barh(ind-0.4, self.data.Procstat['events'][-1,2*i], width, color='r')
                if self.data.Procstat['events'][-1,2*i+1] > 0:
                    plt.barh(ind-0.6, self.data.Procstat['events'][-1,2*i+1], width, color='b')
                if net_freq > 0:
                    plt.barh(ind-0.8, net_freq, width, color='g')
                ylabels.append(self.data.Reactions['Input'][i]['Name'])                
                yvals.append(ind-0.6)                
                ind = ind - 1

        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('frequency',size=24)
        plt.xscale('log')
        plt.yticks(yvals, ylabels)
        plt.show()
    
    def ComputeTOF(self,Product):                       # return TOF and TOF error
        
        # Find the index of the product species
        product_ind = -1       
        for i in enumerate(self.data.Species['gas_spec']):
            if i[1] == Product:
                product_ind = i[0]
        
        # Make sure the index has been found
        if product_ind == -1:
            print 'Product species not found'
        else:
            product_ind = product_ind + self.data.Species['n_surf']         # Adjust index to account for surface species   
        
        
        nRxns = len(self.data.Reactions['Nu'])        
        TOF_contributions = [0 for i in range(nRxns)]              # number of product molecules produced in each reaction        
        for i, elem_stoich in enumerate(self.data.Reactions['Nu']):
            TOF_stoich = elem_stoich[product_ind]
            r = self.data.Binary['propCounter'][-1,i] / self.data.Specnum['t'][-1]      # ergodic average
#            r = self.data.Binary['prop'][-1,i]                                           # non-ergodic average
            TOF_contributions[i] = TOF_stoich * r         
               
        TOF = np.sum(TOF_contributions)
        TOF_fracs = TOF_contributions / TOF             # will need this for sensitivity analysis
#        return TOF
        return {'TOF': TOF, 'TOF_fracs': TOF_fracs}        
        
"""


      
    def BuildEmptyFolders(self):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        if not os.path.isdir(sysinfo['Path']['LocalRunDir']):
            os.mkdir(sysinfo['Path']['LocalRunDir'])
            os.mkdir(sysinfo['Path']['LocalRunDir'] + 'Build/')
            os.mkdir(sysinfo['Path']['LocalRunDir'] + 'Run/')
            os.mkdir(sysinfo['Path']['LocalRunDir'] + 'IntermediateRuns/')
        
        
    def AdjustRuntime(self,WallTime=24*3600):
        MaxLen = MS().MaxdataEntries()
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
        MaxLen = MS().MaxdataEntries()
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
        
    def KMCSeed(self):
        data = random.randint(1,int(2 ** 31 - 1))
        return data
        
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