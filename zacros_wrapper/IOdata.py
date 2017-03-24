import os
import numpy as np
import re
import linecache
import random

from Helper import *
from Lattice import Lattice     # will fill in this object later when it is used
# GRW Comment
class IOdata(object):
    
    '''
    Handles the reading and writing of Zacros input/output files
    '''
    
    def __init__(self):
        
        '''
        Set default values of class variables
        '''
        
        super(IOdata, self).__init__()        
        
        self.Path = ''

# -------------------- Input data --------------------
        
        self.Conditions                       = {}
        self.Conditions['T']                  = None
        self.Conditions['P']                  = None
        self.Conditions['Seed']               = None
        self.Conditions['restart']            = False
        self.Conditions['SimTime']            = {}
        self.Conditions['SimTime']['Max']     = None
        self.Conditions['SimTime']['Actual']  = None
        self.Conditions['WallTime']           = {}
        self.Conditions['WallTime']['Max']    = None
        self.Conditions['WallTime']['Actual'] = None
        self.Conditions['CPUTime']            = None
        self.Conditions['nEvents']            = None
        self.Conditions['MaxStep']            = None
    
        self.Species                          = {}
        self.Species['n_gas']                 = None
        self.Species['gas_spec']              = None
        self.Species['gas_eng']               = None
        self.Species['gas_MW']                = None
        self.Species['gas_molfrac']           = None
        self.Species['n_surf']                = None
        self.Species['surf_spec']             = None
        self.Species['surf_dent']             = None
    
        self.Report                           = {}
        self.Report['specnum']                = [None, None]
        self.Report['procstat']               = [None, None]
        self.Report['hist']                   = [None, None]
        self.Report['event']                  = 'off'
    
        self.Cluster                          = {}
        self.Cluster['nCluster']              = None
        self.Cluster['nClusterVariant']       = None
        self.Cluster['Input']                 = None
        
        self.Reactions                        = {}
        self.Reactions['nrxns']               = None      # does not include reversible reactions for irreversible reactions
        self.Reactions['nrxns_total']         = None      # does include reversible reactions for irreversible reactions
        self.Reactions['Input']               = None
        self.Reactions['names']               = []
        self.Reactions['reversibilities']     = []
        self.Reactions['is_reversible']       = []
        
        self.StateInput                       = {}
        self.StateInput['Type']               = 'none'
        self.StateInput['Struct']             = None
        
        self.Lattice                          = {}
        self.Lattice['Input']                 = None
        self.KMC_lat = Lattice()

        self.scaledown_factors    = None

# -------------------- Output data --------------------

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
        
        self.History                          = []
        self.n_snapshots = 0
        self.n_sites = 0
        self.snap_times = []
        
        self.Binary                           = {}
        self.Binary['cluster']                = ''
        self.Binary['prop']                   = ''
        self.Binary['propCounter']            = ''
        self.Binary['W_sen_anal']             = ''  
        
        self.Reactions['Nu']      = ''
        self.Reactions['UniqNu']  = ''
        
        self.Performance = {}
        self.Performance['t_final']  = 0
        self.Performance['events_occurred']  = 0
        self.Performance['CPU_time']  = 0

# ------------------------------------- Read input files ----------------------------

    def ReadAllInput(self):   
    
        '''
        Read all input files. state_input.dat will be read only if it is present.
        '''
    
        self.ReadSimIn()
        self.ReadLatticeIn()
        self.ReadEngIn()
        self.ReadMechIn()
       
        if os.path.isfile(os.path.join(self.Path, 'state_input.dat')):
            self.ReadStateInput()
    
    def ReadEngIn(self):

        '''
        Read energetics_input.dat
        '''
    
        RawTxt = ReadWithoutBlankLines(os.path.join(self.Path, 'energetics_input.dat'), CommentLines=False)
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
    
        '''
        Read lattice_input.dat
        '''
    
        self.Lattice['Input'] = []
        with open(os.path.join(self.Path, 'lattice_input.dat'),'r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            self.Lattice['Input'].append(i.split('\n')[0])
    
    def ReadStateInput(self):

        '''
        Read state_input.dat
        '''
    
        self.StateInput['Struct'] = []
        with open(os.path.join(self.Path, 'state_input.dat'),'r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            self.StateInput['Struct'].append(i.split('\n')[0])
        self.StateInput['Type'] = 'StateInput'
    
    def ReadMechIn(self):

        '''
        Read mechanism_input.dat
        '''
    
        RawTxt = ReadWithoutBlankLines(os.path.join(self.Path, 'mechanism_input.dat'),CommentLines=True)
        nLines = len(RawTxt)
        StiffCorrLine = -1
        
        nMech = 0
        for i in range(0,nLines):
            if RawTxt[i].split()[0]=='reversible_step' or RawTxt[i].split()[0]=='step':
                nMech += 1
            elif re.search('# Automated stiffness reconditioning employed',RawTxt[i]):
                StiffCorrLine = i
                
        if StiffCorrLine != -1:
            self.scaledown_factors = [np.float(i) for i in RawTxt[StiffCorrLine+2].split(':')[1].split()]
        
        # Initialize varibles
        self.Reactions['nrxns'] = 0
        self.Reactions['nrxns_total'] = 0
        self.Reactions['names'] = []
        self.Reactions['reversibilities'] = []
        self.Reactions['is_reversible'] = []
        
        MechInd = np.array([[0,0]]*nMech)
        Count = 0
        for i in range(0,nLines):
        
            if RawTxt[i].split()[0]=='reversible_step':
                MechInd[Count,0] = i
                self.Reactions['reversibilities'].append(True)
                
            elif RawTxt[i].split()[0]=='step':
                MechInd[Count,0] = i
                self.Reactions['reversibilities'].append(False)
                
            elif RawTxt[i].split()[0]=='end_reversible_step':
                MechInd[Count,1] = i
                Count += 1
                
            elif RawTxt[i].split()[0]=='end_step':
                MechInd[Count,1] = i
                Count += 1

        MechDict = [{'Name':'','nSites':0,'neighboring':'','initial':'','final':'','variant':'','gas_reacs_prods':''} for k in range(0,nMech)]
        for j in range(0,nMech):
        
            reversible = self.Reactions['reversibilities'][j]
        
            # Count the variants
        
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
#                elif not InVariant and i not in StateLine:
#                    print 'Unparsed line in mechanism input:'
#                    print RawTxt[i]
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
            
            
            # Enumerate reversibilities
            
            self.Reactions['nrxns'] += nVariant
            self.Reactions['is_reversible'] += [reversible for i in range(nVariant)]
            
            if reversible:
                self.Reactions['nrxns_total'] += 2 * nVariant
            else:
                self.Reactions['nrxns_total'] += nVariant
            
            
            for k in range(0,nVariant):
                for i in range(variantInd[k,0],variantInd[k,1]):
                    if RawTxt[i].split()[0]=='variant':
                        MechDict[j]['variant'][k]['Name'] = RawTxt[i].split()[1]
                        self.Reactions['names'].append(RawTxt[MechInd[j,0]].split()[1] + '_' + RawTxt[i].split()[1])
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
#                    else:
#                        print 'Unparsed line in mechanism variant:'
#                        print RawTxt[i]
        
        if StiffCorrLine == -1:
            self.scaledown_factors = np.ones(self.Reactions['nrxns'])    
        
        self.Reactions['Input'] = MechDict    
        
    def ReadSimIn(self):
    
        '''
        Read simulation_input.dat
        '''
    
        with open(os.path.join(self.Path, 'simulation_input.dat'),'r') as txt:
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
                        if i.split()[1] == 'infinity':
                            self.Conditions['SimTime']['Max'] = 'inf'
                        else:
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
#                    else:
#                        print 'Unparsed line in simulation_input.dat:'
#                        print i
        
    def StateInc(self,i):
    
        '''
        Read state_input.dat
        '''
    
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

# ------------------------------------- Write input files ---------------------------- '''

    def WriteAllInput(self):
        
        '''
        Write all Zacros input files
        '''
        
        self.WriteSimIn()
        self.WriteMechanism()
        self.WriteEnergetics()
        self.WriteStateIn()
        self.WriteLattice()

    def WriteEnergetics(self):
    
        '''
        Write energetics_input.dat
        '''
    
        nCluster = self.Cluster['nCluster']

        with open(os.path.join(self.Path, 'energetics_input.dat'), 'w') as txt:
            txt.write('energetics\n\n')
            for i in range(0,nCluster):
                txt.write('#'*80 + '\n\n')
                ClusterStr = self.Cluster['Input'][i]
                txt.write('cluster ' + ClusterStr['Name'] + '\n\n')
                txt.write('  sites ' + str(ClusterStr['nSites']) + '\n')
                
                if ClusterStr['neighboring'] != '':
                    txt.write('  neighboring')
                    for j in range(0,len(ClusterStr['neighboring'])):
                        txt.write(' ' + ClusterStr['neighboring'][j])
                    txt.write('\n')
                    
                txt.write('  lattice_state\n')
                for j in range(0,int(ClusterStr['nSites'])):
                    txt.write(ClusterStr['latstate'][j] + '\n')
                    
                nVariant = len(ClusterStr['variant'])              
                txt.write('\n')
                for j in range(0,nVariant):
                    txt.write('  variant ' + ClusterStr['variant'][j]['Name'] + '\n')
                    txt.write(PadStr('    site_types',25))
                    for k in range(0,len(ClusterStr['variant'][j]['site_types'])):
                        txt.write(ClusterStr['variant'][j]['site_types'][k] + ' ')
                    txt.write('\n')
                    if int(ClusterStr['variant'][j]['graph_multiplicity'])>0:
                        txt.write(PadStr('    graph_multiplicity',25) + str(ClusterStr['variant'][j]['graph_multiplicity']) + '\n')
                    txt.write(PadStr('    cluster_eng',25) + str(ClusterStr['variant'][j]['eng']) + '\n')
                    txt.write('  end_variant\n\n')
#                
                txt.write('end_cluster\n\n')
                txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_energetics')
            
    def WriteLattice(self):
    
        '''
        Write lattice_input.dat
        '''
    
        with open(os.path.join(self.Path, 'lattice_input.dat'), 'w') as txt:
            for i in self.Lattice['Input']:
                txt.write(i + '\n')
    
# Need to fix this method
        
    def WriteMechanism(self):
    
        '''
        Write mechanism_input.dat
        '''
    
        if self.scaledown_factors is None:
            SDBool = False
        else:
            SDBool = True
        nMech = len(self.Reactions['Input'])
        StiffCorrCounter = -1
        with open(os.path.join(self.Path, 'mechanism_input.dat'), 'w') as txt:
            txt.write('mechanism\n\n')
            if SDBool:
                txt.write('# Automated stiffness reconditioning employed\n')
                txt.write('# \n')
                txt.write('# SDF: ')
                for i in self.scaledown_factors:
                    txt.write('{0:.5e} \t'.format(i))
                txt.write('\n\n')
            for i in range(0,nMech):
            
                reversible = self.Reactions['reversibilities'][i]
            
                txt.write('#'*80 + '\n\n')
                MechStr = self.Reactions['Input'][i]
                
                if reversible:
                    txt.write('reversible_step ' + MechStr['Name'] + '\n')
                else:
                    txt.write('step ' + MechStr['Name'] + '\n')
                    
                txt.write('  sites ' + str(MechStr['nSites']) + '\n')             
                if MechStr['neighboring'] != '':
                    txt.write('  neighboring')
                    for j in range(0,len(MechStr['neighboring'])):
                        txt.write(' ' + MechStr['neighboring'][j])
                    txt.write('\n')

                if MechStr['gas_reacs_prods'] != '':
                    txt.write('  gas_reacs_prods ' + MechStr['gas_reacs_prods'][0] + ' ' + str(MechStr['gas_reacs_prods'][1]) + '\n')                    
                
                txt.write('  initial\n')
                for j in range(0,int(MechStr['nSites'])):
                    txt.write(MechStr['initial'][j] + '\n')
                    
                txt.write('  final\n')
                for j in range(0,int(MechStr['nSites'])):
                    txt.write(MechStr['final'][j] + '\n')
                    
                nVariant = len(MechStr['variant'])
                txt.write('\n')
                for j in range(nVariant):
                    txt.write('  variant ' + MechStr['variant'][j]['Name'] + '\n')
                    txt.write(PadStr('    site_types',25))
                    for k in range(0,len(MechStr['variant'][j]['site_types'])):
                        txt.write(MechStr['variant'][j]['site_types'][k] + ' ')
                    txt.write('\n')
                    pre_exp = MechStr['variant'][j]['pre_expon']
                    if SDBool:
                        StiffCorrCounter += 1
                        if self.scaledown_factors[StiffCorrCounter] != 1:                            
                            txt.write(PadStr('    pre_expon',25) 
                            + '{0:.5e}'.format(pre_exp))
                            txt.write('    # Pre-exponential has been rescaled by a factor of ' + 
                            '{0:.5e}'.format(self.scaledown_factors[StiffCorrCounter]) + ' \n')
                        else:
                            txt.write(PadStr('    pre_expon',25) 
                            + '{0:.5e}'.format(pre_exp) + '\n')
                    else:                        
                        txt.write(PadStr('    pre_expon',25) 
                        + '{0:.5e}'.format(pre_exp) + '\n')
                    if reversible:
                        txt.write(PadStr('    pe_ratio',25) + '{0:.5e}'.format(MechStr['variant'][j]['pe_ratio']) + '\n')
                    txt.write(PadStr('    activ_eng',25) + str(MechStr['variant'][j]['activ_eng']) + '\n')
                    if MechStr['variant'][j]['prox_factor'] != '':
                        txt.write(PadStr('    prox_factor',25) + str(MechStr['variant'][j]['prox_factor']) + '\n')
                    txt.write('  end_variant\n\n')
                
                if reversible:
                    txt.write('end_reversible_step\n\n')
                else:
                    txt.write('end_step\n\n')
                    
                txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_mechanism')
    
    def WriteSimIn(self):
    
        '''
        Write simulation_input.dat
        '''
    
        with open(os.path.join(self.Path, 'simulation_input.dat'), 'w') as txt:
            SeedTxt = ''
            if self.Conditions['Seed'] == '':
                self.Conditions['Seed'] = random.randint(10000, 99999)
                SeedTxt = '      #Random seed from Python wrapper'
            
            txt.write('#KMC simulation specification\n\n')
            txt.write('random_seed' + ' '*9 + str(self.Conditions['Seed']) + SeedTxt + '\n\n')
            
            txt.write('temperature         ' + '{0:.5f}'.format(self.Conditions['T']) + '\n')
            txt.write('pressure            ' + '{0:.5f}'.format(self.Conditions['P']) + '\n\n') 
            txt.write('n_gas_species       ' + str(self.Species['n_gas']) + '\n')
            txt.write('gas_specs_names     ')
            for i in range(0,self.Species['n_gas']):
                txt.write(PadStr(self.Species['gas_spec'][i],14) + ' ')
                
            GasList  = ['gas_energies','gas_molec_weights','gas_molar_fracs']
            GasList2 = ['gas_eng','gas_MW','gas_molfrac']
            for j in range(0,len(GasList)):
                txt.write('\n' + PadStr(GasList[j],19) + ' ')
                for i in range(0,self.Species['n_gas']):
                    txt.write(PadStr( '{0:.5f}'.format(self.Species[GasList2[j]][i]),14) + ' ')

            txt.write('\n\n')
            txt.write('n_surf_species      ' + str(self.Species['n_surf']) + '\n')
            txt.write('surf_specs_names    ')
            for i in range(0,self.Species['n_surf']):
                txt.write(PadStr(self.Species['surf_spec'][i],14) + ' ')
            txt.write('\nsurf_specs_dent     ')
            for i in range(0,self.Species['n_surf']):
                txt.write(PadStr(str(self.Species['surf_dent'][i]),15))
            txt.write('\n\n')    
            
            if self.Report['hist'][0] == 'off':
                txt.write('snapshots           off\n')
            elif self.Report['hist'][0] == 'event':
                txt.write('snapshots           on '  + self.Report['hist'][0] + ' ' + str(int(self.Report['hist'][1])) + '\n')
            elif self.Report['hist'][0] == 'time':
                txt.write('snapshots           on '  + self.Report['hist'][0] + ' ' + str(np.float(self.Report['hist'][1])) + '\n')
                
            if self.Report['procstat'][0] == 'off':
                txt.write('process_statistics  off\n')
            elif self.Report['procstat'][0] == 'event':
                txt.write('process_statistics  on '  + self.Report['procstat'][0] + ' ' + str(int(self.Report['procstat'][1])) + '\n')
            elif self.Report['procstat'][0] == 'time':
                txt.write('process_statistics  on '  + self.Report['procstat'][0] + ' ' + str(np.float(self.Report['procstat'][1])) + '\n')
                
            if self.Report['specnum'][0] == 'off':
                txt.write('species_numbers     off\n')
            elif self.Report['specnum'][0] == 'event':
                txt.write('species_numbers     on '  + self.Report['specnum'][0] + ' ' + str(int(self.Report['specnum'][1])) + '\n')
            elif self.Report['specnum'][0] == 'time':
                txt.write('species_numbers     on '  + self.Report['specnum'][0] + ' ' + str(np.float(self.Report['specnum'][1])) + '\n')

            txt.write('event_report ' + self.Report['event'] + '\n\n')
            if self.Conditions['MaxStep'] == '' or re.search('inf',str(self.Conditions['MaxStep'])):
                txt.write('max_steps           infinity\n')
            else:
                txt.write('max_steps           ' + str(self.Conditions['MaxStep']) + '\n')
            
            if self.Conditions['SimTime']['Max'] == '' or re.search('inf',str(self.Conditions['SimTime']['Max'])):
                txt.write('max_time            infinity\n')
            else:
                txt.write('max_time            ' + str(self.Conditions['SimTime']['Max']) + '\n')
            if self.Conditions['WallTime']['Max'] == '' or re.search('inf',str(self.Conditions['WallTime']['Max'])):
                txt.write('\n')
#                txt.write('\nwall_time           ' + str(3600*24*365*10) + '\n\n')      # 10 years
            else:
                txt.write('\nwall_time           ' + str(self.Conditions['WallTime']['Max']) + '\n\n')
            
            if not self.Conditions['restart']:
                txt.write('no_restart\n')
            txt.write('\nfinish\n')

    def WriteStateIn(self):
    
        '''
        Write state_input.dat
        '''
    
        if self.StateInput['Type'] == 'none':
            return
            
        if self.StateInput['Type'] == 'StateInput':   #Copy from prior state_input file
            with open(os.path.join(self.Path, 'state_input.dat'), 'w') as txt:
                for i in self.StateInput['Struct']:
                    txt.write(i + '\n')
                    
        elif self.StateInput['Type'] == 'history':   #Copy from prior history_output file
        
            Lattice = self.StateInput['Struct']
            UniqSpec = np.unique(Lattice[np.not_equal(Lattice[:,2],0),1])
            nAds = len(UniqSpec)
            SpecIden = [0] * nAds
            AdsInfo = [[] for i in range(0,nAds)]
            DentInfo = [[] for i in range(0,nAds)]
            for i in range(0,nAds):
                for j in range(0,Lattice.shape[0]):
                    if UniqSpec[i] == Lattice[j,1]:
                        AdsInfo[i].append(j+1)
                        DentInfo[i].append(Lattice[j,3])
                        SpecIden[i] = Lattice[j,2]
            
            if nAds > 0:
                with open(os.path.join(self.Path, 'state_input.dat'),'w') as txt:
                    txt.write('initial_state\n')
                    for i in range(0,nAds):
                        txt.write('  seed_on_sites  ' + PadStr(self.Species['surf_spec'][SpecIden[i]-1],10))
                        for j in range(0,len(DentInfo[i])):
                            for k in range(0,len(DentInfo[i])):
                                if j + 1 == DentInfo[i][k]:
                                    txt.write(str(AdsInfo[i][k]) + '  ')
                        txt.write('\n')
                    txt.write('end_initial_state\n')
        else:
            print 'Unrecognized state_input type'
            print 'state_input not written'
                
#------------------------------------- Read output files ----------------------------

    def ReadAllOutput(self, build_lattice = False):
        
        '''
        Read all Zacros output files
        Set build_lattice = True if you want to build the lattice object file.
        This will make it take a lot longer to read
        '''
        
        self.ReadAllInput()
        
        if self.CheckComplete():
            
            # Standard output files
            self.ReadGeneral()
            self.ReadProcstat()
            self.ReadSpecnum()
            if build_lattice:
                self.KMC_lat.Read_lattice_output(os.path.join(self.Path, 'lattice_output.txt'))
            self.ReadHistory()
            
            # Extra binary files
            if os.path.isfile(os.path.join(self.Path, 'PropCounter_output.bin')):            
                self.ReadProp(1)            
            if os.path.isfile(os.path.join(self.Path, 'SA_output.bin')):
                self.ReadSA()           
#            self.ReadCluster()
#            self.ReadProp(0)            
            
        else:
            print 'general_output.txt not found in ' + self.Path
  
  
    def CheckComplete(self):
    
        '''
        Check to see if a Zacros run has completed successfully
        '''
    
        Complete = False
        if os.path.isfile(os.path.join(self.Path, 'general_output.txt')):
            with open(os.path.join(self.Path, 'general_output.txt'),'r') as txt:
                RawTxt = txt.readlines()
            for i in RawTxt:
                if re.search('Normal termination',i):
                    Complete = True
        return Complete

    def ReadGeneral(self):
    
        '''
        Read general_output.txt
        '''
    
        with open(os.path.join(self.Path, 'general_output.txt'),'r') as txt:
            RawTxt = txt.readlines()                
                
        for i in range(0,len(RawTxt)):
            if re.search('Number of elementary steps:',RawTxt[i]):
                nRxn = np.int(RawTxt[i].split(':')[1])
            elif re.search('Current KMC time:',RawTxt[i]):
                self.Performance['t_final'] = np.float(RawTxt[i].split(':')[1])
            elif re.search('Events occurred:',RawTxt[i]):
                self.Performance['events_occurred'] = np.float(RawTxt[i].split(':')[1])
            elif re.search('Elapsed CPU time:',RawTxt[i]):
                after_colon = RawTxt[i].split(':')[1]
                self.Performance['CPU_time'] = np.float(after_colon.split(' ')[-2])
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
            nu = [0] * (self.Species['n_surf'] + self.Species['n_gas'])
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
                        SurfInd = [k for k in range(0,len(self.Species['surf_spec'])) if SurfIden == self.Species['surf_spec'][k]][0]
                        nu[SurfInd] += Sign
                elif RxnStrList[j] != '->' and RxnStrList[j] != '+':
                    GasInd = [k for k in range(0,len(self.Species['gas_spec'])) if RxnStrList[j] == self.Species['gas_spec'][k]][0]
                    nu[self.Species['n_surf'] + GasInd] += Sign
            nuList.append(nu)

        self.Reactions['Nu']      = nuList
        self.Reactions['UniqNu']  = ReturnUnique(nuList).tolist()    
    
    def ReadHistory(self):
        
        '''
        Read history_output.txt
        '''
        
        # Check if file exists
        if not os.path.isfile(os.path.join(self.Path, 'history_output.txt')):
            return
        
        # Count number of sites
        with open(os.path.join(self.Path, 'lattice_output.txt'),'r') as txt:
            RawTxt = txt.readlines()
        nSites = len(RawTxt) - 2
        self.n_sites = nSites
        HistPath = os.path.join(self.Path, 'history_output.txt')
        nLines = rawbigcount(HistPath)
        self.n_snapshots = (nLines-6)/(nSites+2)
        self.History = []
        self.snap_times = []
        
        for snap_ind in range(self.n_snapshots):
            snap_data = np.array([[0]*4]*nSites)
            linecache.clearcache()
            snap_header = linecache.getline(HistPath, 8 + snap_ind * (nSites+2)-1).split()
            self.snap_times.append(snap_header[3])
            for i in range(0,nSites):               
                snap_data[i,:] = linecache.getline(HistPath, 8 + snap_ind * (nSites+2)+i).split()
            self.History.append(snap_data)
        
    def ReadProcstat(self):
    
        '''
        Read procstat_output.txt
        '''
    
        MaxLen = np.int(2e4)
        with open(os.path.join(self.Path, 'procstat_output.txt'),'r') as txt:
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

        '''
        Read specnum_output.txt
        '''
    
        MaxLen = np.int(2e4)
        with open(os.path.join(self.Path, 'specnum_output.txt'),'r') as txt:
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
    
        '''
        Read clusterocc.bin
        '''
    
        dt=np.dtype(np.int32)
        virtual_arr = np.memmap(os.path.join(self.Path, 'clusterocc.bin'), dt, "r")
        nCluster = self.Cluster['nClusterVariant']
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nCluster)
        virtual_arr = virtual_arr[:nNum]
        self.Binary['cluster'] = np.array(np.reshape(virtual_arr,[nNum/nCluster,nCluster])[::self.Specnum['Spacing']])
        del virtual_arr
    
    def ReadProp(self,Mode):
    
        '''
        Mode = 0: Read Prop_output.bin
        Mode = 1: Read PropCounter_output.bin
        '''
    
        dt=np.dtype(np.float64)
        if Mode==0:     #Instantaneous propensities
            FileName = 'Prop_output.bin'
            
        elif Mode==1:   #Integral propensities
            FileName = 'PropCounter_output.bin'
        
        virtual_arr = np.memmap(os.path.join(self.Path, FileName), dt, "r")
        nRxn = len(self.Reactions['Nu'])
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
    
        '''
        Read SA_output.bin
        '''
    
        dt=np.dtype(np.float64)
        FileName = 'SA_output.bin'
        if os.path.isfile(os.path.join(self.Path, FileName)):
            virtual_arr = np.memmap(os.path.join(self.Path, FileName), dt, "r")
            nRxn = len(self.Reactions['Nu'])
            nNum = virtual_arr.shape[0]
            nNum = nNum - (nNum % nRxn)
            virtual_arr = virtual_arr[:nNum]
            self.Binary['W_sen_anal'] = np.reshape(virtual_arr,[nNum/nRxn,nRxn])
            self.Binary['W_sen_anal'] = np.array(self.Binary['W_sen_anal'][::self.Specnum['Spacing']])
            
            del virtual_arr
        else:
            print 'No sensitivity analysis output file'