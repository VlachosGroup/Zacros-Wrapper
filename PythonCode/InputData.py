# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 14:37:02 2016

@author: robieta
"""

import GeneralUtilities as ut
import KMCUtilities as KMCut
import numpy as np
import os
import re

class :
    def __init__(self):
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
    
    def ReadAllInput(self,Path,Cnd = KMCut()):
        
        Cnd = self.ReadSimIn(Path,Cnd)
        Cnd = self.ReadLatticeIn(Path,Cnd)
        Cnd = self.ReadEngIn(Path,Cnd)
        Cnd = self.ReadMechIn(Path,Cnd)
       
        if os.path.isfile(Path + 'state_input.dat'):
            Cnd = self.ReadStateInput(Path,Cnd)
        return Cnd
      
    def ReadEngIn(self,Path,Cnd = KMCut()): 
        RawTxt = ut.GeneralUtilities().ReadWithoutBlankLines(Path + 'energetics_input.dat',CommentLines=False)
        nLines = len(RawTxt)
        
        nCluster = 0
        for i in range(0,nLines):
            if RawTxt[i].split()[0]=='cluster':
                nCluster += 1
                
        ClusterInd = np.array([[0,0]]*nCluster)
        Count = 0
        for i in range(0,nLines):
            if RawTxt[i].split()[0]=='cluster':
                ClusterInd[Count,0] = i
            if RawTxt[i].split()[0]=='end_cluster':
                ClusterInd[Count,1] = i
                Count += 1
        ClusterDict = [{'Name':'','nSites':0,'neighboring':'','latstate':'','variant':''} for k in range(0,nCluster)]
        
        nClusterTotal = 0
        for j in range(0,nCluster):
            ClusterDict[j]['Name'] = RawTxt[ClusterInd[j,0]].split()[1]
            Count = 0
            for i in range(ClusterInd[j,0]+1,ClusterInd[j,1]):
                if RawTxt[i].split()[0]=='variant':
                    Count += 1
                elif RawTxt[i].split()[0]=='sites':
                    nSites = int(RawTxt[i].split()[1])
                    ClusterDict[j]['nSites'] = nSites
                elif RawTxt[i].split()[0]=='neighboring':
                    neighbor = RawTxt[i].split()[1:]
                    ClusterDict[j]['neighboring']=neighbor
                elif RawTxt[i].split()[0]=='lattice_state':
                    LatState = RawTxt[i+1:i+1+nSites]
                    ClusterDict[j]['latstate']=LatState
                    for k in range(0,len(ClusterDict[j]['latstate'])):
                            ClusterDict[j]['latstate'][k] = ClusterDict[j]['latstate'][k].split('\n')[0]
                    
            nVariant = Count
            nClusterTotal += nVariant
            ClusterDict[j]['variant']=[{'Name':'','site_types':'','graph_multiplicity':0,'eng':0.} for k in range(0,nVariant)]
            variantInd = np.array([[0,0]]*nVariant)
            Count = 0
            for i in range(ClusterInd[j,0]+1,ClusterInd[j,1]):
                if RawTxt[i].split()[0]=='variant':
                    variantInd[Count,0] = i
                if RawTxt[i].split()[0]=='end_variant':
                    variantInd[Count,1] = i
                    Count +=1
                    
            for k in range(0,nVariant):
                for i in range(variantInd[k,0],variantInd[k,1]):
                    if RawTxt[i].split()[0]=='variant':
                        ClusterDict[j]['variant'][k]['Name'] = RawTxt[i].split()[1]
                    elif RawTxt[i].split()[0]=='site_types':
                        ClusterDict[j]['variant'][k]['site_types'] = RawTxt[i].split()[1:]
                    elif RawTxt[i].split()[0]=='graph_multiplicity':
                        ClusterDict[j]['variant'][k]['graph_multiplicity'] = RawTxt[i].split()[1]
                    elif RawTxt[i].split()[0]=='cluster_eng':
                        ClusterDict[j]['variant'][k]['eng'] = RawTxt[i].split()[1]
        
        Cnd['Cluster']['Input'] = ClusterDict
        Cnd['Cluster']['nCluster'] = len(ClusterDict)
        Cnd['Cluster']['nClusterVariant'] = nClusterTotal
        return Cnd
        
    def ReadLatticeIn(self,Path,Cnd = KMCut()):
        Cnd['Lattice']['Input'] = []
        with open(Path + 'lattice_input.dat','r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            Cnd['Lattice']['Input'].append(i.split('\n')[0])
        return Cnd
    
    def ReadStateInput(self,Path,Cnd = KMCut()): 
        Cnd['StateInput']['Struct'] = []
        with open(Path + 'state_input.dat','r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            Cnd['StateInput']['Struct'].append(i.split('\n')[0])
        Cnd['StateInput']['Type'] = 'StateInput'
        
        return Cnd
    
    def ReadMechIn(self,Path,Cnd = KMCut()): 
        RawTxt = ut.GeneralUtilities().ReadWithoutBlankLines(Path + 'mechanism_input.dat',CommentLines=True)
        nLines = len(RawTxt)
        StiffCorrLine = -1
        
        nMech = 0
        for i in range(0,nLines):
            if RawTxt[i].split()[0]=='reversible_step':
                nMech += 1
            elif RawTxt[i].split()[0]=='step':
                raise NameError('Wrapper does not support irreversable steps')
            elif re.search('# Automated stiffness reconditioning employed',RawTxt[i]):
                StiffCorrLine = i
                
        if StiffCorrLine != -1:
            Cnd['StiffnessRecondition']['Mode'] = RawTxt[StiffCorrLine+1].split(':')[1].split('\n')[0].split()[0]
            Cnd['StiffnessRecondition']['APSdF'] = [np.float(i) for i in RawTxt[StiffCorrLine+2].split(':')[1].split()]
        
        MechInd = np.array([[0,0]]*nMech)
        Count = 0
        for i in range(0,nLines):
            if RawTxt[i].split()[0]=='reversible_step':
                MechInd[Count,0] = i
            if RawTxt[i].split()[0]=='end_reversible_step':
                MechInd[Count,1] = i
                Count += 1

        MechDict = [{'Name':'','nSites':0,'neighboring':'','initial':'','final':'','variant':'','gas_reacs_prods':''} for k in range(0,nMech)]
        for j in range(0,nMech):
            MechDict[j]['Name'] = RawTxt[MechInd[j,0]].split()[1]
            Count = 0
            InVariant = False
            StateLine = []
            for i in range(MechInd[j,0]+1,MechInd[j,1]):
                if RawTxt[i].split()[0]=='variant':
                    Count += 1   
                    InVariant = True
                elif RawTxt[i].split()[0]=='end_variant':
                    InVariant = False
                elif RawTxt[i].split()[0]=='gas_reacs_prods':    
                    MechDict[j]['gas_reacs_prods'] = RawTxt[i].split()[1:]
                elif RawTxt[i].split()[0]=='sites':
                    nSites = int(RawTxt[i].split()[1])
                    MechDict[j]['nSites'] = nSites
                elif RawTxt[i].split()[0]=='neighboring':
                    neighbor = RawTxt[i].split()[1:]
                    MechDict[j]['neighboring']=neighbor
                elif RawTxt[i].split()[0]=='initial':
                    LatState = RawTxt[i+1:i+1+nSites]
                    MechDict[j]['initial']=LatState
                    for k in range(0,len(MechDict[j]['initial'])):
                            MechDict[j]['initial'][k] = MechDict[j]['initial'][k].split('\n')[0]
                    for k in range(0,nSites):
                        StateLine.append(i+1+k)
                elif RawTxt[i].split()[0]=='final':
                    LatState = RawTxt[i+1:i+1+nSites]
                    MechDict[j]['final']=LatState
                    for k in range(0,len(MechDict[j]['initial'])):
                            MechDict[j]['final'][k] = MechDict[j]['final'][k].split('\n')[0]
                    for k in range(0,nSites):
                        StateLine.append(i+1+k)
                elif not InVariant and i not in StateLine:
                    print 'Unparsed line in mechanism input:'
                    print RawTxt[i]
            nVariant = Count  
            MechDict[j]['variant']=[{'Name':'','site_types':'','pre_expon':'','pe_ratio':'','activ_eng':'','prox_factor':''} for k in range(0,nVariant)]     
            variantInd = np.array([[0,0]]*nVariant)
            Count = 0
            for i in range(MechInd[j,0]+1,MechInd[j,1]):
                if RawTxt[i].split()[0]=='variant':
                    variantInd[Count,0] = i
                if RawTxt[i].split()[0]=='end_variant':
                    variantInd[Count,1] = i
                    Count +=1
            for k in range(0,nVariant):
                for i in range(variantInd[k,0],variantInd[k,1]):
                    if RawTxt[i].split()[0]=='variant':
                        MechDict[j]['variant'][k]['Name'] = RawTxt[i].split()[1]
                    elif RawTxt[i].split()[0]=='site_types':
                        MechDict[j]['variant'][k]['site_types'] = RawTxt[i].split()[1:]
                    elif RawTxt[i].split()[0]=='pre_expon':
                        MechDict[j]['variant'][k]['pre_expon'] = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0]=='pe_ratio':
                        MechDict[j]['variant'][k]['pe_ratio'] = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0]=='activ_eng':
                        MechDict[j]['variant'][k]['activ_eng'] = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0]=='prox_factor':
                        MechDict[j]['variant'][k]['prox_factor'] = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == '#':
                        pass
                    else:
                        print 'Unparsed line in mechanism variant:'
                        print RawTxt[i]
        
        Cnd['Reactions']['Input'] = MechDict
        return Cnd        
        
    def ReadSimIn(self,Path,Cnd = KMCut()):
        with open(Path + 'simulation_input.dat','r') as txt:
            RawTxt = txt.readlines()
            
        Cnd['Conditions']['restart'] = True
        for i in RawTxt:
            if len(i.split())>0:
                if i[0] != '#':
                    i=i.split('#')[0] # Don't parse comments
                    if i.split()[0] == 'temperature':
                        Cnd['Conditions']['T']          = np.float(i.split()[1])
                    elif i.split()[0] == 'pressure':
                        Cnd['Conditions']['P']          = np.float(i.split()[1])
                    elif i.split()[0] == 'random_seed':
                        Cnd['Conditions']['Seed']       = np.int(i.split()[1])
                    elif i.split()[0] == 'no_restart':
                        Cnd['Conditions']['restart'] = False
                    elif i.split()[0] == 'gas_specs_names':
                        Cnd['Species']['gas_spec']      = i.split()[1:]
                        Cnd['Species']['n_gas']         = len(Cnd['Species']['gas_spec'])
                    elif i.split()[0] == 'gas_energies':
                        Cnd['Species']['gas_eng']       = []
                        for j in i.split()[1:]:
                            Cnd['Species']['gas_eng'].append(np.float(j))
                    elif i.split()[0] == 'gas_molec_weights':
                        Cnd['Species']['gas_MW']        = []
                        for j in i.split()[1:]:
                            Cnd['Species']['gas_MW'].append(np.float(j))
                    elif i.split()[0] == 'gas_molar_fracs':
                        Cnd['Species']['gas_molfrac']   = []
                        for j in i.split()[1:]:
                            Cnd['Species']['gas_molfrac'].append(np.float(j))
                    elif i.split()[0] == 'surf_specs_names':
                        Cnd['Species']['surf_spec']     = i.split()[1:]
                        Cnd['Species']['n_surf']        = len(Cnd['Species']['surf_spec'])
                    elif i.split()[0] == 'surf_specs_dent':
                        Cnd['Species']['surf_dent']     = []
                        for j in i.split()[1:]:
                            Cnd['Species']['surf_dent'].append(np.int(j))
                    
                    elif i.split()[0] == 'event_report':
                        Cnd['Report']['event'] = i.split()[1]
                    elif i.split()[0] == 'snapshots':
                        Cnd['Report']['hist']           = self.StateInc(i)
                    elif i.split()[0] == 'process_statistics':
                        Cnd['Report']['procstat']       = self.StateInc(i) 
                    elif i.split()[0] == 'species_numbers':
                        Cnd['Report']['specnum']        = self.StateInc(i)
                    
                    elif i.split()[0] == 'max_time':
                        Cnd['Conditions']['SimTime']['Max'] = np.float(i.split()[1])
                    elif i.split()[0] == 'max_steps':
                        if i.split()[1] == 'infinity':
                            Cnd['Conditions']['MaxStep'] = 'inf'
                        else:
                            Cnd['Conditions']['MaxStep'] = int(i.split()[1])
                    elif i.split()[0] == 'wall_time':
                        Cnd['Conditions']['WallTime']['Max'] = np.int(i.split()[1])
                    elif i.split()[0] == 'finish' or i.split()[0] == 'n_gas_species' or i.split()[0] == 'n_surf_species':
                        pass
                    else:
                        print 'Unparsed line in simulation_input.dat:'
                        print i
            
        return Cnd
        
    def StateInc(self,i):
        if re.search('off',i):
            state = 'off'
            inc = ''
        elif re.search('on time',i):
            state = 'time'
            inc = np.float(i.split()[3])
        elif re.search('on event',i):
            state = 'event'
            inc = np.int(i.split()[3])
        Output = [state,inc]
        return Output