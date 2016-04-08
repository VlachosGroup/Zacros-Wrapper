# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 13:32:13 2016

@author: robieta
"""

import GeneralUtilities as ut
import KMCUtilities as KMCut
import linecache
import numpy as np
import os
import pickle
import re
import ReadInputFiles as RI
import sys

sys.path.append("..")
from DynamicFiles.MachineSpecifics import MachineSpecifics as MS

class ReadOutputFiles:
    def __init__(self):
        pass
    
    def ReadJobOutput(self,Path):
        DirList = ut.GeneralUtilities().GetDir(Path)

        nDir = len(DirList)
        if os.path.isfile(Path + 'CndList.p'):
            CndList = pickle.load(open( Path + 'CndList.p', "rb" ))
        else:
            CndList = ['' for i in range(nDir)]
            for i in range(nDir):
                RunPath = Path + DirList[i] + '/'
                CndList[i]  = KMCut.KMCUtilities().InitializeCnd()
                CndList[i] = RI.ReadInputFiles().ReadAll(RunPath,CndList[i])
                CndList[i] = self.ReadAll(RunPath,CndList[i])
                print str(i+1) + ' / ' + str(nDir)
            pickle.dump( CndList, open( Path + 'CndList.p', "wb" ) )
        return CndList
            
        
            
    
    def ReadAll(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        # Note on computational efficiency:
        #    Testing indicates that python passes arguements to functions
        #    using pointers. (expected) This means that a large Cnd dictionary
        #    should not significantly impede performance even if it is passed
        #    between functions multiple times.
    
        Cnd = self.GetGeneral(Path,Cnd)
        Cnd = self.GetProcstat(Path,Cnd)
        Cnd = self.GetSpecnum(Path,Cnd)
        Cnd = self.GetBinaries(Path,Cnd)
        Cnd = self.GetHistory(Path,Cnd)
        Cnd = KMCut.KMCUtilities().CndToTau(Cnd)
        return Cnd
    
    def GetGeneral(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        if os.path.isfile(Path + 'general_output.txt'):
            Cnd = self.ReadGeneral(Path,Cnd)
        return Cnd
        
    def GetHistory(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        if os.path.isfile(Path + 'history_output.txt'):
            Cnd = self.ReadHistory(Path,Cnd)
        return Cnd
    
    def GetProcstat(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        if os.path.isfile(Path + 'procstat_output.txt'):
            Cnd = self.ReadProcstat(Path,Cnd)
        return Cnd
        
    def GetSpecnum(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        if os.path.isfile(Path + 'specnum_output.txt'):
            Cnd = self.ReadSpecnum(Path,Cnd)
        return Cnd
        
    def GetBinaries(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        if Cnd['Cluster'] != '' and Cnd['Reactions'] != '' and Cnd['Specnum'] != '' and Cnd['Procstat'] != '':
            if os.path.isfile(Path + 'clusterocc.bin'):
                self.ReadCluster(Path,Cnd)
            if os.path.isfile(Path + 'Prop_output.bin'):
                self.ReadProp(Path,Cnd,0)
            if os.path.isfile(Path + 'PropCounter_output.bin'):
                self.ReadProp(Path,Cnd,1)
        return Cnd
  
    def ReadGeneral(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        with open(Path + 'general_output.txt','r') as txt:
            RawTxt = txt.readlines()

        Cnd['Conditions']['restart'] = 'true'
        for i in RawTxt:
            if re.search('Temperature:',i):
                Cnd['Conditions']['T']          = np.float(i.split(':')[1])
            elif re.search('Pressure:',i):
                Cnd['Conditions']['P']          = np.float(i.split(':')[1])
            elif re.search('Random sequence with seed:',i):
                Cnd['Conditions']['Seed']       = np.int(i.split(':')[1])
            elif re.search('Keyword no_restart parsed.',i):
                Cnd['Conditions']['restart'] = 'false'
            

            elif re.search('Gas species names:',i):
                Cnd['Species']['gas_spec']      = i.split(':')[1].split()
                Cnd['Species']['n_gas']         = len(Cnd['Species']['gas_spec'])
            elif re.search('Gas species energies:',i):
                Cnd['Species']['gas_eng']       = []
                for j in i.split(':')[1].split():
                    Cnd['Species']['gas_eng'].append(np.float(j))
            elif re.search('Gas species molecular weights:',i):
                Cnd['Species']['gas_MW']        = []
                for j in i.split(':')[1].split():
                    Cnd['Species']['gas_MW'].append(np.float(j))
            elif re.search('Gas species molar fractions:',i):
                Cnd['Species']['gas_molfrac']   = []
                for j in i.split(':')[1].split():
                    Cnd['Species']['gas_molfrac'].append(np.float(j))
                    
            elif re.search('Surface species names:',i):
                Cnd['Species']['surf_spec']     = i.split(':')[1].split()
                Cnd['Species']['n_surf']        = len(Cnd['Species']['surf_spec'])
            elif re.search('Surface species dentation:',i):
                Cnd['Species']['surf_dent']     = []
                for j in i.split(':')[1].split():
                    Cnd['Species']['surf_dent'].append(np.int(j))
            
            elif re.search('Snapshots will be saved',i) or re.search('Snapshot saving',i):
                Cnd['Report']['hist']           = self.StateInc(i)
            elif re.search('Process statistics will be reported',i) or re.search('Process statistics reporting',i):
                Cnd['Report']['procstat']       = self.StateInc(i)         
            elif re.search('Species number will be reported',i) or re.search('Species numbers reporting',i):
                Cnd['Report']['specnum']        = self.StateInc(i)
            elif re.search('Event reporting',i):
                Cnd['Report']['event']          = i.split()[-1]
                
            elif re.search('Max simulated time:',i):
                if re.search('maximum allowed value',i):
                    Cnd['Conditions']['SimTime']['Max'] = 'infinity'
                else:
                    Cnd['Conditions']['SimTime']['Max'] = np.float(i.split()[-1])
            elif re.search('Current KMC time:',i):
                Cnd['Conditions']['SimTime']['Actual'] = np.float(i.split(':')[1])
            elif re.search('Allowed walltime in seconds:',i):
                Cnd['Conditions']['WallTime']['Max'] = np.int(i.split(':')[1])
            elif re.search('Elapsed clock time:',i):
                Cnd['Conditions']['WallTime']['Actual'] = np.float(i.split(':')[1].split()[0])
            elif re.search('Elapsed CPU time:',i):
                Cnd['Conditions']['CPUTime'] = np.float(i.split(':')[1].split()[0])
            elif re.search('Events occurred:',i):
                Cnd['Conditions']['nEvents'] = np.int(i.split(':')[1])
                
            elif re.search('Number of clusters:',i):
                Cnd['Cluster']['nClusterVariant'] = np.int(i.split(':')[1])
                
                
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
            nu = [0] * (Cnd['Species']['n_surf'] + Cnd['Species']['n_gas'])
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
                        SurfInd = [k for k in range(0,len(Cnd['Species']['surf_spec'])) if SurfIden == Cnd['Species']['surf_spec'][k]][0]
                        nu[SurfInd] += Sign
                elif RxnStrList[j] != '->' and RxnStrList[j] != '+':
                    GasInd = [k for k in range(0,len(Cnd['Species']['gas_spec'])) if RxnStrList[j] == Cnd['Species']['gas_spec'][k]][0]
                    nu[Cnd['Species']['n_surf'] + GasInd] += Sign
            nuList.append(nu)

        Cnd['Reactions']['Names']   = RxnNameList   
        Cnd['Reactions']['Nu']      = nuList
        Cnd['Reactions']['UniqNu']  = ut.GeneralUtilities().ReturnUnique(nuList).tolist()
                
        return Cnd
 
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
    
    def ReadHistory(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        Cnd['History']['Final'] = self.ReadSnapshot(Path,-1)           
        return Cnd
             
    def ReadSnapshot(self,Path,Snapshot):
        with open(Path + 'lattice_output.txt','r') as txt:
            RawTxt = txt.readlines()
        nSites = len(RawTxt) - 2
        SnapshotArray = np.array([[0]*4]*nSites)
        HistPath = Path + 'history_output.txt'
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
     
        
    def ReadProcstat(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        MaxLen = MS().MaxOutputEntries()    
        with open(Path + 'procstat_output.txt','r') as txt:
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
        
        Cnd['Procstat']['Spacing'] = Spacing
        Cnd['Procstat']['t'] = t
        Cnd['Procstat']['events'] = events
        
        return Cnd
    
    def ReadSpecnum(self,Path,Cnd = KMCut.KMCUtilities().InitializeCnd()):
        MaxLen = MS().MaxOutputEntries()
        with open(Path + 'specnum_output.txt','r') as txt:
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
        
        Cnd['Specnum']['Spacing'] = Spacing 
        Cnd['Specnum']['nEvents']   = nEvents
        Cnd['Specnum']['t']         = t
        Cnd['Specnum']['T']         = T
        Cnd['Specnum']['E']         = E
        Cnd['Specnum']['spec']      = spec

        return Cnd
    
    def ReadCluster(self,Path,Cnd):
        dt=np.dtype(np.int32)
        virtual_arr = np.memmap(Path + 'clusterocc.bin', dt, "r")
        nCluster = Cnd['Cluster']['nClusterVariant']
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nCluster)
        virtual_arr = virtual_arr[:nNum]
        Cnd['Binary']['cluster'] = np.array(np.reshape(virtual_arr,[nNum/nCluster,nCluster])[::Cnd['Specnum']['Spacing']])
        del virtual_arr
        return Cnd
    
    def ReadProp(self,Path,Cnd,Mode):
        dt=np.dtype(np.float64)
        if Mode==0:     #Instantaneous propensities
            FileName = 'Prop_output.bin'
            
        elif Mode==1:   #Integral propensities
            FileName = 'PropCounter_output.bin'
        
        virtual_arr = np.memmap(Path + FileName, dt, "r")
        nRxn = len(Cnd['Reactions']['Nu'])
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]
            
        if Mode==0:
            Cnd['Binary']['prop'] = np.reshape(virtual_arr,[nNum/nRxn,nRxn])
            Cnd['Binary']['prop'] = np.array(Cnd['Binary']['prop'][::Cnd['Procstat']['Spacing']])
        if Mode==1:
            Cnd['Binary']['propCounter'] = np.reshape(virtual_arr,[nNum/nRxn,nRxn])
            Cnd['Binary']['propCounter'] = np.array(Cnd['Binary']['propCounter'][::Cnd['Procstat']['Spacing']])
        
        del virtual_arr
        return Cnd
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    