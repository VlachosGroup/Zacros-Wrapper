# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 13:32:13 2016

@author: robieta
"""

import os
import re
import numpy as np
import GeneralUtilities as ut
import linecache
from InputData import InputData

class OutputData:
    
    def __init__(self):

        self.Path = ''

        self.input = InputData()

        self.Specnum                          = {}
        self.Specnum['Spacing']               = ''
        self.Specnum['t']                     = ''
        self.Specnum['nEvents']               = ''
        self.Specnum['T']                     = ''
        self.Specnum['E']                     = ''
        self.Specnum['spec']                  = ''
        
        self.Procstat                         = {}
        self.Procstat['Spacing']              = ''
        self.Procstat['t']                    = ''
        self.Procstat['events']               = ''
        
        self.History                          = {}
        self.History['Final']                 = ''
        
        self.Binary                           = {}
        self.Binary['cluster']                = ''
        self.Binary['prop']                   = ''
        self.Binary['propCounter']            = ''
        self.Binary['W_sen_anal']             = ''  
          
    def ReadAllOutput(self):
        
        self.input.Path = self.Path
        self.input.ReadAllInput()
        
        print 'Reading output files in ' + self.Path
        if self.CheckComplete():
            self.ReadGeneral()
            self.ReadProcstat()
            self.ReadSpecnum()
#            self.ReadHistory()
            
            self.ReadCluster()
            self.ReadProp(0)            
            self.ReadProp(1)            
            self.ReadSA()
        else:
            print 'general_output.txt not found in ' + self.Path
  
    def CheckComplete(self):
        Complete = False
        if os.path.isfile(self.Path + 'general_output.txt'):
            with open(self.Path + 'general_output.txt','r') as txt:
                RawTxt = txt.readlines()
            for i in RawTxt:
                if re.search('Normal termination',i):
                    Complete = True
        return Complete

    def ReadGeneral(self):
        with open(self.Path + 'general_output.txt','r') as txt:
            RawTxt = txt.readlines()                
                
        for i in range(0,len(RawTxt)):
            if re.search('Number of elementary steps:',RawTxt[i]):
                nRxn = np.int(RawTxt[i].split(':')[1])
            elif re.search('Reaction network:',RawTxt[i]):
                RxnStartLine = i + 2
                
        if RawTxt[RxnStartLine].split()[0] == '1.':
            NameInd = 1
        else:
            NameInd = 0
        
        RxnNameList = []
        nuList = []
        for i in range(RxnStartLine,RxnStartLine + nRxn):
            RxnName = RawTxt[i].split()[NameInd][:-1]
            RxnNameList.append(RxnName)
            RxnStr = RawTxt[i][re.search('Reaction:', RawTxt[i]).end():]
            RxnStrList = RxnStr.split()
            nu = [0] * (self.input.Species['n_surf'] + self.input.Species['n_gas'])
            for j in range(0,len(RxnStrList)):
                if RxnStrList[j] == '->':
                    ArrowInd = j
            for j in range(0,len(RxnStrList)):
                if j < ArrowInd:
                    Sign = -1
                else:
                    Sign = 1
                    
                if re.search('\(',RxnStrList[j]):
                    SurfIden = re.sub(r'\([^)]*\)', '', RxnStrList[j])                    
                    if SurfIden != '*':
                        SurfInd = [k for k in range(0,len(self.input.Species['surf_spec'])) if SurfIden == self.input.Species['surf_spec'][k]][0]
                        nu[SurfInd] += Sign
                elif RxnStrList[j] != '->' and RxnStrList[j] != '+':
                    GasInd = [k for k in range(0,len(self.input.Species['gas_spec'])) if RxnStrList[j] == self.input.Species['gas_spec'][k]][0]
                    nu[self.input.Species['n_surf'] + GasInd] += Sign
            nuList.append(nu)

        self.input.Reactions['Names']   = RxnNameList
        self.input.Reactions['Nu']      = nuList
        self.input.Reactions['UniqNu']  = ut.GeneralUtilities().ReturnUnique(nuList).tolist()
                

    def StateInc(self,i):
        if re.search('turned off',i):
            state = 'off'
            inc = ''
        elif re.search('time units',i):
            state = 'time'
            inc = np.float(i.split()[-3])
        elif re.search('events',i):
            state = 'event'
            inc = np.int(i.split()[-2])
        Output = [state,inc]
        return Output
    
    def ReadHistory(self):
        self.History['Final'] = self.ReadSnapshot(-1)           
             
    def ReadSnapshot(self,Snapshot):
        with open(self.Path + 'lattice_output.txt','r') as txt:
            RawTxt = txt.readlines()
        nSites = len(RawTxt) - 2
        SnapshotArray = np.array([[0]*4]*nSites)
        HistPath = self.Path + 'history_output.txt'
        nLines = ut.GeneralUtilities().rawbigcount(HistPath)
        nSnapshot = np.float(nLines-6)/(nSites+2)
        if nSnapshot != int(nSnapshot):
            raise ValueError('Index error in the history_state.txt read')
        if Snapshot < 0:
            Snapshot = int(nSnapshot) + Snapshot
        linecache.clearcache()
        for i in range(0,nSites):
            SnapshotArray[i,:] = linecache.getline(HistPath, 8+Snapshot*(nSites+2)+i).split()
        return SnapshotArray
     
        
    def ReadProcstat(self):
        MaxLen = np.int(2e4)
        with open(self.Path + 'procstat_output.txt','r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen * 3: # Procstat uses 3 lines per outputs
            Spacing = np.int(np.floor((len(RawTxt)-1)/(MaxLen*3)))
            RawTxt2 = []
            for i in range(0,MaxLen):
                RawTxt2.append(RawTxt[i*Spacing*3+1])
                RawTxt2.append(RawTxt[i*Spacing*3+2])
                RawTxt2.append(RawTxt[i*Spacing*3+3])
        else:
            Spacing = 1
            RawTxt2 = RawTxt[1:]
            
            
        t = []
        events = []
        for i in range(0,len(RawTxt2)/3):
            t.append(np.float(RawTxt2[i*3].split()[3]))
            eventsTemp = RawTxt2[i*3+2].split()[1:]
            for j in range(0,len(eventsTemp)):
                eventsTemp[j] = np.int(eventsTemp[j])
            events.append(eventsTemp)
        
        self.Procstat['Spacing'] = Spacing
        self.Procstat['t'] = np.asarray(t)
        self.Procstat['events'] = np.asarray(events)
    
    def ReadSpecnum(self):
        MaxLen = np.int(2e4)
        with open(self.Path + 'specnum_output.txt','r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen:
            Spacing = np.int(np.floor((len(RawTxt)-1)/MaxLen))
            RawTxt2 = []
            for i in range(0,MaxLen):
                RawTxt2.append(RawTxt[i*Spacing+1])
        else:
            Spacing = 1
            RawTxt2 = RawTxt[1:]
        
        nEvents = []
        t = []
        T = []
        E = [] 
        spec = []
        
        for i in range(0,len(RawTxt2)):
            LineSplit = RawTxt2[i].split()
            nEvents.append(np.int(LineSplit[1]))
            t.append(np.float(LineSplit[2]))
            T.append(np.float(LineSplit[3]))
            E.append(np.float(LineSplit[4]))
            specTemp = LineSplit[5:]
            for j in range(0,len(specTemp)):
                specTemp[j] = np.int(specTemp[j])
            spec.append(specTemp)
        
        self.Specnum['Spacing'] = Spacing 
        self.Specnum['nEvents']   = np.asarray(nEvents)
        self.Specnum['t']         = np.asarray(t)
        self.Specnum['T']         = np.asarray(T)
        self.Specnum['E']         = np.asarray(E)
        self.Specnum['spec']      = np.asarray(spec)
    
    def ReadCluster(self):
        dt=np.dtype(np.int32)
        virtual_arr = np.memmap(self.Path + 'clusterocc.bin', dt, "r")
        nCluster = self.input.Cluster['nClusterVariant']
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nCluster)
        virtual_arr = virtual_arr[:nNum]
        self.Binary['cluster'] = np.array(np.reshape(virtual_arr,[nNum/nCluster,nCluster])[::self.Specnum['Spacing']])
        del virtual_arr
    
    def ReadProp(self,Mode):
        dt=np.dtype(np.float64)
        if Mode==0:     #Instantaneous propensities
            FileName = 'Prop_output.bin'
            
        elif Mode==1:   #Integral propensities
            FileName = 'PropCounter_output.bin'
        
        virtual_arr = np.memmap(self.Path + FileName, dt, "r")
        nRxn = len(self.input.Reactions['Nu'])
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]
            
        if Mode==0:
            self.Binary['prop'] = np.reshape(virtual_arr,[nNum/nRxn,nRxn])
            self.Binary['prop'] = np.array(self.Binary['prop'][::self.Procstat['Spacing']])
        if Mode==1:
            self.Binary['propCounter'] = np.reshape(virtual_arr,[nNum/nRxn,nRxn])
            self.Binary['propCounter'] = np.array(self.Binary['propCounter'][::self.Procstat['Spacing']])
        
        del virtual_arr
    
    def ReadSA(self):
        dt=np.dtype(np.float64)
        FileName = 'SA_output.bin'
        if os.path.isfile(self.Path + FileName):
            virtual_arr = np.memmap(self.Path + FileName, dt, "r")
            nRxn = len(self.input.Reactions['Nu'])
            nNum = virtual_arr.shape[0]
            nNum = nNum - (nNum % nRxn)
            virtual_arr = virtual_arr[:nNum]
            self.Binary['W_sen_anal'] = np.reshape(virtual_arr,[nNum/nRxn,nRxn])
            self.Binary['W_sen_anal'] = np.array(self.Binary['W_sen_anal'][::self.Specnum['Spacing']])  
            del virtual_arr
        else:
            print 'No sensitivity analysis output file'