# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 14:37:02 2016

@author: robieta
"""

import os
import GeneralUtilities as ut
import numpy as np
import re

class InputData:
    
    def __init__(self):
        
        self.Path                             = ''
        
        self.Conditions                       = {}
        self.Conditions['T']                  = ''
        self.Conditions['P']                  = ''
        self.Conditions['Seed']               = ''
        self.Conditions['restart']            = ''
        self.Conditions['SimTime']            = {}
        self.Conditions['SimTime']['Max']     = ''
        self.Conditions['SimTime']['Actual']  = ''
        self.Conditions['WallTime']           = {}
        self.Conditions['WallTime']['Max']    = ''
        self.Conditions['WallTime']['Actual'] = ''
        self.Conditions['CPUTime']            = ''
        self.Conditions['nEvents']            = ''
        self.Conditions['MaxStep']            = ''
    
        self.Species                          = {}
        self.Species['n_gas']                 = ''
        self.Species['gas_spec']              = ''
        self.Species['gas_eng']               = ''
        self.Species['gas_MW']                = ''
        self.Species['gas_molfrac']           = ''
        self.Species['n_surf']                = ''
        self.Species['surf_spec']             = ''
        self.Species['surf_dent']             = ''
    
        self.Report                           = {}
        self.Report['specnum']                = ['','']
        self.Report['procstat']               = ['','']
        self.Report['hist']                   = ['','']
        self.Report['event']                  = ''
    
        self.Cluster                          = {}
        self.Cluster['nCluster']              = ''
        self.Cluster['nClusterVariant']       = ''
        self.Cluster['Input']                 = ''
        
        self.Reactions                        = {}
        self.Reactions['nrxns']               = ''
        self.Reactions['Input']               = ''
        self.Reactions['Names']               = ''
        self.Reactions['Nu']                  = ''
        self.Reactions['UniqNu']              = ''        
        
        self.StateInput                       = {}
        self.StateInput['Type']               = ''
        self.StateInput['Struct']             = ''
        
        self.Lattice                          = {}
        self.Lattice['Input']                 = ''

        # Will deal with this later
        self.StiffnessRecondition             = {}
        self.StiffnessRecondition['Mode']     = ''
        self.StiffnessRecondition['APSdF']    = ''

    def ReadAllInput(self):
    
        print 'Reading input files in ' + self.Path    
    
        self.ReadSimIn()
        self.ReadLatticeIn()
        self.ReadEngIn()
        self.ReadMechIn()
       
        if os.path.isfile(self.Path + 'state_input.dat'):
            self.ReadStateInput()
    
    def ReadEngIn(self): 
        RawTxt = ut.GeneralUtilities().ReadWithoutBlankLines(self.Path + 'energetics_input.dat',CommentLines=False)
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
        
        self.Cluster['Input'] = ClusterDict
        self.Cluster['nCluster'] = len(ClusterDict)
        self.Cluster['nClusterVariant'] = nClusterTotal
        
    def ReadLatticeIn(self):
        self.Lattice['Input'] = []
        with open(self.Path + 'lattice_input.dat','r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            self.Lattice['Input'].append(i.split('\n')[0])
    
    def ReadStateInput(self): 
        self.StateInput['Struct'] = []
        with open(self.Path + 'state_input.dat','r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            self.StateInput['Struct'].append(i.split('\n')[0])
        self.StateInput['Type'] = 'StateInput'
    
    def ReadMechIn(self): 
        RawTxt = ut.GeneralUtilities().ReadWithoutBlankLines(self.Path + 'mechanism_input.dat',CommentLines=True)
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
            self.StiffnessRecondition['Mode'] = RawTxt[StiffCorrLine+1].split(':')[1].split('\n')[0].split()[0]
            self.StiffnessRecondition['APSdF'] = [np.float(i) for i in RawTxt[StiffCorrLine+2].split(':')[1].split()]
        
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
        
        self.Reactions['Input'] = MechDict    
        
        # Find unique reactions
#        Cnd['Reactions']['Names']   = RxnNameList   
#        Cnd['Reactions']['Nu']      = nuList
#        Cnd['Reactions']['UniqNu']  = ut.GeneralUtilities().ReturnUnique(nuList).tolist()
        
    def ReadSimIn(self):
        with open(self.Path + 'simulation_input.dat','r') as txt:
            RawTxt = txt.readlines()
            
        self.Conditions['restart'] = True
        for i in RawTxt:
            if len(i.split())>0:
                if i[0] != '#':
                    i=i.split('#')[0] # Don't parse comments
                    if i.split()[0] == 'temperature':
                        self.Conditions['T']          = np.float(i.split()[1])
                    elif i.split()[0] == 'pressure':
                        self.Conditions['P']          = np.float(i.split()[1])
                    elif i.split()[0] == 'random_seed':
                        self.Conditions['Seed']       = np.int(i.split()[1])
                    elif i.split()[0] == 'no_restart':
                        self.Conditions['restart'] = False
                    elif i.split()[0] == 'gas_specs_names':
                        self.Species['gas_spec']      = i.split()[1:]
                        self.Species['n_gas']         = len(self.Species['gas_spec'])
                    elif i.split()[0] == 'gas_energies':
                        self.Species['gas_eng']       = []
                        for j in i.split()[1:]:
                            self.Species['gas_eng'].append(np.float(j))
                    elif i.split()[0] == 'gas_molec_weights':
                        self.Species['gas_MW']        = []
                        for j in i.split()[1:]:
                            self.Species['gas_MW'].append(np.float(j))
                    elif i.split()[0] == 'gas_molar_fracs':
                        self.Species['gas_molfrac']   = []
                        for j in i.split()[1:]:
                            self.Species['gas_molfrac'].append(np.float(j))
                    elif i.split()[0] == 'surf_specs_names':
                        self.Species['surf_spec']     = i.split()[1:]
                        self.Species['n_surf']        = len(self.Species['surf_spec'])
                    elif i.split()[0] == 'surf_specs_dent':
                        self.Species['surf_dent']     = []
                        for j in i.split()[1:]:
                            self.Species['surf_dent'].append(np.int(j))
                    
                    elif i.split()[0] == 'event_report':
                        self.Report['event'] = i.split()[1]
                    elif i.split()[0] == 'snapshots':
                        self.Report['hist']           = self.StateInc(i)
                    elif i.split()[0] == 'process_statistics':
                        self.Report['procstat']       = self.StateInc(i) 
                    elif i.split()[0] == 'species_numbers':
                        self.Report['specnum']        = self.StateInc(i)
                    
                    elif i.split()[0] == 'max_time':
                        self.Conditions['SimTime']['Max'] = np.float(i.split()[1])
                    elif i.split()[0] == 'max_steps':
                        if i.split()[1] == 'infinity':
                            self.Conditions['MaxStep'] = 'inf'
                        else:
                            self.Conditions['MaxStep'] = int(i.split()[1])
                    elif i.split()[0] == 'wall_time':
                        self.Conditions['WallTime']['Max'] = np.int(i.split()[1])
                    elif i.split()[0] == 'finish' or i.split()[0] == 'n_gas_species' or i.split()[0] == 'n_surf_species':
                        pass
                    else:
                        print 'Unparsed line in simulation_input.dat:'
                        print i
        
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
