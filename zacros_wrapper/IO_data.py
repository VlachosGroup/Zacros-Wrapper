
import os
import numpy as np
import re as _re
import random as _random
import linecache as _linecache

from utils import constant as _c
from utils import ReadWithoutBlankLines as _ReadWithoutBlankLines
from utils import ReturnUnique as _ReturnUnique
from utils import rawbigcount as _rawbigcount

from Thermochemistry.io_.excel import read_excel
from Thermochemistry.models.empirical.zacros import Zacros

'''
 Class definitions for:

     Zacros input file      Class name
     -----------------      ------------
     simulation_input.dat   SimIn
     energetics_input.dat   ClusterIn
     mechanism_input.dat    MechanismIn
     
     lattice_input.dat      LatticeIn
     state_input.dat        StateIn

     Zacros output file    Class name
     ------------------    ------------
     general_output.txt    Performance
     history_output.txt    History
     procstat_output.txt   Procstat
     specnum_output.txt    SpecnumOut
'''


'''
============ Classes to handle input files ============
'''


class SimIn():

    '''
    Handles input from simulation_input.dat
    '''

    fname = 'simulation_input.dat'
    
    def __init__(self):
    
        self.TPD = False            # flag for temperature programmed desorption (TPD) mode
        self.TPD_start = None
        self.TPD_ramp = None
        self.T = None
        self.P = None
        
        self.gas_spec = []
        self.n_gas = 0
        self.gas_eng = []
        self.gas_MW = []
        self.gas_molfrac = []
        self.surf_spec = []
        self.surf_dent = []
        
        self.Seed = None
        self.restart = True
        self.WallTime_Max = None
        self.MaxStep = None
        self.SimTime_Max = None
        self.hist = None
        
        
    def ReadIn(self, fldr):
    
        '''
        Read simulation_input.dat
        '''
        
        with open(os.path.join(fldr, self.fname), 'r') as txt:
            RawTxt = txt.readlines()

        for i in RawTxt:
            if len(i.split()) > 0:
                if i[0] != '#':
                    i = i.split('#')[0]  # Don't parse comments
                    if i.split()[0] == 'temperature':
                        if i.split()[1] == 'ramp':
                            self.TPD = True
                            self.TPD_start = np.float(i.split()[2])
                            self.TPD_ramp = np.float(i.split()[3])
                        else:
                            self.T = np.float(i.split()[1])
                    elif i.split()[0] == 'pressure':
                        self.P = np.float(i.split()[1])
                    elif i.split()[0] == 'random_seed':
                        self.Seed = np.int(i.split()[1])
                    elif i.split()[0] == 'no_restart':
                        self.restart = False
                    elif i.split()[0] == 'gas_specs_names':
                        self.gas_spec = i.split()[1:]
                        self.n_gas = len(self.gas_spec)
                    elif i.split()[0] == 'gas_energies':
                        for j in i.split()[1:]:
                            self.gas_eng.append(np.float(j))
                    elif i.split()[0] == 'gas_molec_weights':
                        for j in i.split()[1:]:
                            self.gas_MW.append(np.float(j))
                    elif i.split()[0] == 'gas_molar_fracs':
                        for j in i.split()[1:]:
                            self.gas_molfrac.append(np.float(j))
                    elif i.split()[0] == 'surf_specs_names':
                        self.surf_spec = i.split()[1:]
                        self.n_surf = len(self.surf_spec)
                    elif i.split()[0] == 'surf_specs_dent':
                        for j in i.split()[1:]:
                            self.surf_dent.append(np.int(j))
                    elif i.split()[0] == 'event_report':
                        self.event = i.split()[1]
                    elif i.split()[0] == 'snapshots':
                        self.hist = StateInc(i)
                    elif i.split()[0] == 'process_statistics':
                        self.procstat = StateInc(i)
                    elif i.split()[0] == 'species_numbers':
                        self.specnum = StateInc(i)
                    elif i.split()[0] == 'max_time':
                        if i.split()[1] == 'infinity':
                            self.SimTime_Max = 'inf'
                        else:
                            self.SimTime_Max =\
                             np.float(i.split()[1])
                    elif i.split()[0] == 'max_steps':
                        if i.split()[1] == 'infinity':
                            self.MaxStep = 'inf'
                        else:
                            self.MaxStep = int(i.split()[1])
                    elif i.split()[0] == 'wall_time':
                        self.WallTime_Max = np.int(i.split()[1])
                    elif i.split()[0] == 'finish' or\
                                         i.split()[0] == 'n_gas_species' or\
                                         i.split()[0] == 'n_surf_species':
                        pass
                        

    def WriteIn(self, fldr):
    
        '''
        Write simulation_input.dat
        '''

        with open(os.path.join(fldr, self.fname), 'w') as txt:
            SeedTxt = ''
            if self.Seed is None:
                _random.seed()
                self.Seed = _random.randint(10000, 99999)
                SeedTxt = '# Random seed from Python wrapper'

            txt.write('#KMC simulation specification\n\n')
            txt.write('{:20}{:15}{}\n\n'.format('random_seed',
                      str(self.Seed), SeedTxt))
            
            # Write out temperature, which depends on TPD or constant temperature mode
            if self.TPD:
                txt.write('{:20}{:5.1f}{:5.1f}\n'.format('temperature\t ramp', self.TPD_start, self.TPD_ramp))
            else:
                txt.write('{:20}{:5.1f}\n'.format('temperature', self.T))

            txt.write('{:20}{:5.5e}\n\n'.format('pressure', self.P))
            txt.write('{:20}{}\n'.format('n_gas_species',
                      str(self.n_gas)))
            txt.write('{:20}'.format('gas_specs_names'))
            for i in range(0, self.n_gas):
                txt.write('{:15}'.format(self.gas_spec[i]))

            GasList = ['gas_energies', 'gas_molec_weights', 'gas_molar_fracs']
            GasList2 = ['gas_eng', 'gas_MW', 'gas_molfrac']
            for j in range(0, len(GasList)):
                txt.write('\n{:20}'.format(GasList[j]))
                for i in range(0, self.n_gas):
                    txt.write('{:15}'.format(str(getattr(self,
                              GasList2[j])[i])))

            txt.write('\n\n{:20}{}\n'.format('n_surf_species',
                      str(self.n_surf)))
            txt.write('{:20}'.format('surf_specs_names'))
            for i in range(0, self.n_surf):
                txt.write('{:15}'.format(self.surf_spec[i]))
            txt.write('\n{:20}'.format('surf_specs_dent'))

            for i in range(0, self.n_surf):
                txt.write('{:15}'.format(str(self.surf_dent[i])))
            txt.write('\n\n')

            if self.hist is None:
                txt.write('{:20}{}\n'.format('snapshots', 'off'))
            elif self.hist[0] == 'event':
                txt.write('{:20}{} {} {}\n'.format('snapshots', 'on',
                          self.hist[0],
                          str(int(self.hist[1]))))
            elif self.hist[0] == 'time':
                txt.write('{:20}{} {} {}\n'.format('snapshots', 'on',
                          self.hist[0],
                          str(np.float(self.hist[1]))))
                          
            if self.procstat[0] == 'off':
                txt.write('process_statistics  off\n')
            elif self.procstat[0] == 'event':
                txt.write('{:20}{} {} {}\n'.format('process_statistics', 'on',
                          self.procstat[0],
                          str(int(self.procstat[1]))))
            elif self.procstat[0] == 'time':
                txt.write('{:20}{} {} {}\n'.format('process_statistics', 'on',
                          self.procstat[0],
                          str(np.float(self.procstat[1]))))

            if self.specnum[0] == 'off':
                txt.write('species_numbers     off\n')
            elif self.specnum[0] == 'event':
                txt.write('{:20}{} {} {}\n'.format('species_numbers', 'on',
                          self.specnum[0],
                          str(int(self.specnum[1]))))
            elif self.specnum[0] == 'time':
                txt.write('{:20}{} {} {}\n'.format('species_numbers', 'on',
                          self.specnum[0],
                          str(np.float(self.specnum[1]))))
            txt.write('{:20}{}\n\n'.format('event_report',
                      self.event))

            if self.MaxStep is None or\
               _re.search('inf', str(self.MaxStep)):
                txt.write('{:20}{}\n'.format('max_steps', 'infinity'))
            else:
                txt.write('{:20}{}\n'.format('max_steps',
                          str(self.MaxStep)))

            if self.SimTime_Max is None or\
               _re.search('inf', str(self.SimTime_Max)):
                txt.write('{:20}{}\n'.format('max_time', 'infinity\n'))
            else:
                txt.write('{:20}{}\n'.format('max_time',
                          str(self.SimTime_Max)))
            if self.WallTime_Max is None or\
               _re.search('inf', str(self.WallTime_Max)):
                txt.write('\n')
            else:
                txt.write('\n{:20}{}\n\n'.format('wall_time',
                          str(self.WallTime_Max)))

            if not self.restart:
                txt.write('no_restart\n')
            txt.write('finish\n')
            


class Cluster():

    '''
    Cluster in energetics_input.dat
    '''

    def __init__(self):
    
        self.name = None
        self.variant_list = []
        self.sites = None
        self.neighboring = None
        self.latstate = None
        
    
class cluster_variant():

    '''
    Variant of a cluster in energetics_input.dat
    '''

    def __init__(self):
    
        self.name = 'var'
        self.site_types = None
        self.graph_multiplicity = 1
        self.cluster_eng = 0.0
    

class ClusterIn(object):

    '''
    Handles data from energetics_input.dat
    '''

    fname = 'energetics_input.dat'

    
    def __init__(self):
    
        self.cluster_list = []
        

    def FindCluster(self, Cluster_Num):     # FIX THIS METHOD
    
        '''
        Method finds the Cluster and Variant index of the nth
        Cluster-Variant where n is specified by Cluster_Num and
        the indices are returned as C.index (Cluster)
        and V.index (Variant) such that Cluster[C_index].variant_name[V_index]
        represents the name of the nth Cluster-Variant
        '''
        
        Cluster_Num = int(Cluster_Num)
        Tvariants = sum(s.nVariant for s in self.cluster_list)
        if Tvariants >= Cluster_Num and Cluster_Num >= 1:
            var = []
            for s in self.cluster_list:
                var.append(s.nVariant)
            var = np.array(var)
            C_index = np.argmin(var.cumsum() < Cluster_Num)
            V_index = Cluster_Num - sum(var[0:C_index]) - 1
        else:
            C_index = V_index = -1
        return(C_index, V_index)
        
    
    def get_num_clusters(self):
    
        n_clusters = 0
        for clustr in self.cluster_list:
            for varnt in clustr.variant_list:
                n_clusters += 1
                
        return n_clusters
    
    def ReadIn(self, fldr):
    
        '''
        Read energetics_input.dat
        '''
        
        RawTxt = _ReadWithoutBlankLines(os.path.join(fldr, self.fname), CommentLines=False)
        nLines = len(RawTxt)
    
        nClusters = 0
        for i in range(0, nLines):
            if RawTxt[i].split()[0] == 'cluster':
                nClusters += 1
    
        ClusterInd = np.array([ [0, 0] ] * nClusters)
        Count = 0
        for i in range(0, nLines):
            if RawTxt[i].split()[0] == 'cluster':
                ClusterInd[Count, 0] = i
            if RawTxt[i].split()[0] == 'end_cluster':
                ClusterInd[Count, 1] = i
                Count += 1
    
        nClusterTotal = 0
        self.cluster_list = [Cluster() for j in range(nClusters)]
        
        # Loop through all clusters
        for j in range(nClusters):
        
            self.cluster_list[j].name = RawTxt[ClusterInd[j, 0]].split()[1]
            n_variants = 0
            
            for i in range(ClusterInd[j, 0] + 1, ClusterInd[j, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    n_variants += 1
                elif RawTxt[i].split()[0] == 'sites':
                    self.cluster_list[j].sites = int(RawTxt[i].split()[1])
                elif RawTxt[i].split()[0] == 'neighboring':
                    self.cluster_list[j].neighboring = RawTxt[i].split()[1:]
                elif RawTxt[i].split()[0] == 'lattice_state':
                    self.cluster_list[j].latstate = RawTxt[i + 1:i + 1 +
                                                    self.cluster_list[j].sites]
                    for k in range(0, len(self.cluster_list[j].latstate)):
                        self.cluster_list[j].latstate[k] =\
                        self.cluster_list[j].latstate[k].split('\n')[0]
    
            nClusterTotal += n_variants
            
            # Find beginning and ending lines for each variant
            
            if n_variants == 0:
            
                n_variants = 1
                variantInd = np.array([[0, 0]]*n_variants)
                variantInd[0, 0] = ClusterInd[j, 0] + 1
                variantInd[0, 1] = ClusterInd[j, 1]
                
            else:
            
                variantInd = np.array([[0, 0]]*n_variants)
                Count = 0
                for i in range(ClusterInd[j, 0]+1, ClusterInd[j, 1]):
                    if RawTxt[i].split()[0] == 'variant':
                        variantInd[Count, 0] = i
                    if RawTxt[i].split()[0] == 'end_variant':
                        variantInd[Count, 1] = i
                        Count += 1
                    
            self.cluster_list[j].variant_list = [cluster_variant() for k in range(n_variants)]

            # Loop through all variants for this cluster
            for k in range(n_variants):
            
                for i in range(variantInd[k, 0], variantInd[k, 1]):
                
                    if RawTxt[i].split()[0] == 'variant':
                        self.cluster_list[j].variant_list[k].name = RawTxt[i].split()[1]
                    elif RawTxt[i].split()[0] == 'site_types':
                        self.cluster_list[j].variant_list[k].site_types = RawTxt[i].split()[1:]
                    elif RawTxt[i].split()[0] == 'graph_multiplicity':
                        self.cluster_list[j].variant_list[k].graph_multiplicity = int(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'cluster_eng':
                        self.cluster_list[j].variant_list[k].cluster_eng = float(RawTxt[i].split()[1])
     
                        
    def WriteIn(self, fldr):
    
        '''
        Write energetics_input.dat
        '''
    
        with open(os.path.join(fldr, self.fname), 'w') as txt:
        
            txt.write('energetics\n\n')
            
            for clustr in self.cluster_list:

                txt.write('#'*80 + '\n\n')
                txt.write('cluster ' + clustr.name + '\n\n')
                txt.write('  sites ' + str(clustr.sites) + '\n')
    
                if not clustr.neighboring is None:    # None when there is a point cluster
                    txt.write('  neighboring')
                    for j in range(0, len(clustr.neighboring)):
                        txt.write(' ' + clustr.neighboring[j])
                    txt.write('\n')
    
                txt.write('  lattice_state\n')
                for j in range(0, int(clustr.sites)):
                    txt.write(clustr.latstate[j] + '\n')
    
                txt.write('\n')
                
                for varnt in clustr.variant_list:
                
                    txt.write('  {} {}\n'.format('variant',
                            varnt.name ))
                    txt.write('    {:25}'.format('site_types'))
                    for k in range( len(varnt.site_types )):
                        txt.write('{} '.format
                                (varnt.site_types[k]))
                    txt.write('\n')
                    if int(varnt.graph_multiplicity ) > 0:
                        txt.write('    {:25}{}\n'.format('graph_multiplicity',
                                str(varnt.graph_multiplicity )))
                    txt.write('    {:25}{}\n'.format('cluster_eng',
                            str(varnt.cluster_eng)))
                    txt.write('  end_variant\n\n')
    
                txt.write('end_cluster\n\n')
                
            txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_energetics')
            

class Reaction():

    '''
    Reaction in mechanism_input.dat
    '''

    def __init__(self):
    
        self.name = None
        self.is_reversible = True
        self.gas_reacs_prods = None
        self.sites = None
        self.neighboring = None
        self.initial = None
        self.final = None
        self.variant_list = []

class rxn_variant():

    '''
    Variant of a reaction in mechanism_input.dat
    '''

    def __init__(self):
    
        self.name = 'var'
        self.site_types = None              # site types
        self.pre_expon = None               # pre-exponential factor
        self.pe_ratio = None                # partial equilibrium ratio
        self.activ_eng = 0.0               # activation energy
        self.prox_factor = 0.5
        
        self.scaledown_factor = 1.0
        
            
class MechanismIn(object):

    '''
    Handles input from mechanism_input.dat
    '''
    
    fname = 'mechanism_input.dat'

    def __init__(self):
    
        self.rxn_list = []
        self.include_scaledown = False
        
    
    def get_num_rxns(self):
        
        n_rxns = 0
        for i in self.rxn_list:
            n_rxns += len( i.variant_list )
                
        return n_rxns
        
    
    def get_rxn_var_inds(self, rxn_ind):
    
        ind = 0
        for i in range(len(self.rxn_list)):
            for j in range(len( self.rxn_list[i].variant_list )):
                if ind == rxn_ind:
                    return [i,j]
                ind +=1
                
    
    def FindReaction(self, Reaction_Num):       # FIX THIS METHOD
    
        '''
        Method finds the Reaction and Variant index of the nth
        Reaction-Variant where n is specified by Cluster_Num and
        the indices are returned as R.index (Reaction)
        and V.index (Variant) such that Reaction[R_index].variant_name[V_index]
        represents the name of the nth Reaction-Variant
        '''
        
        Reaction_Num = int(Reaction_Num)
        Tvariants = sum(s.nVariant for s in self.Reaction)
        if Tvariants >= Reaction_Num and Reaction_Num >= 1:
            var = []
            for s in self.Reaction:
                var.append(s.nVariant)
            var = np.array(var)
            R_index = np.argmin(var.cumsum() < Reaction_Num)
            V_index = Reaction_Num - sum(var[0:R_index]) - 1
        else:
            R_index = V_index = -1
        return(R_index, V_index)
        
        
    def ReadIn(self, fldr):
    
        '''
        Read mechanism_input.dat
        '''
        
        RawTxt = _ReadWithoutBlankLines(os.path.join(fldr, self.fname), CommentLines=True)
        nLines = len(RawTxt)
        StiffCorrLine = -1
    
        self.rxn_list = []
        n_rxns = 0
        for i in range(nLines):
            if RawTxt[i].split()[0] == 'reversible_step':
                self.rxn_list.append( Reaction() )
                n_rxns += 1
            elif RawTxt[i].split()[0] == 'step':
                self.rxn_list.append( Reaction() )
                self.rxn_list[-1].is_reversible = False
                n_rxns += 1
            elif _re.search('# Automated stiffness reconditioning employed',
                            RawTxt[i]):
                StiffCorrLine = i
                self.include_scaledown = True
    
        if StiffCorrLine != -1:
            scaledown_factor_list = [np.float(i) for i in RawTxt[StiffCorrLine+2].split(':')[1].split()]
    
        # Identify which lines of text are for each reaction
        MechInd = np.array([[0, 0]]*n_rxns)
        Count = 0
        for i in range(nLines):
    
            if RawTxt[i].split()[0] == 'reversible_step' or RawTxt[i].split()[0] == 'step':
                MechInd[Count, 0] = i
    
            elif RawTxt[i].split()[0] == 'end_reversible_step':
                MechInd[Count, 1] = i
                Count += 1
    
            elif RawTxt[i].split()[0] == 'end_step':
                MechInd[Count, 1] = i
                Count += 1
        
        all_rxn_ind = 0     # Use this index to assign scaledown factors
        
        # Loop over list of recations
        for j in range(n_rxns):
    
            # Count the variants
    
            self.rxn_list[j].name = RawTxt[MechInd[j, 0]].split()[1]
            n_variants = 0
            InVariant = False
            StateLine = []
            for i in range(MechInd[j, 0] + 1, MechInd[j, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    n_variants += 1
                    InVariant = True
                elif RawTxt[i].split()[0] == 'end_variant':
                    InVariant = False
                elif RawTxt[i].split()[0] == 'gas_reacs_prods':
                    self.rxn_list[j].gas_reacs_prods = RawTxt[i].split()[1:]
                elif RawTxt[i].split()[0] == 'sites':
                    nSites = int(RawTxt[i].split()[1])
                    self.rxn_list[j].sites = nSites
                elif RawTxt[i].split()[0] == 'neighboring':
                    self.rxn_list[j].neighboring = RawTxt[i].split()[1:]
                elif RawTxt[i].split()[0] == 'initial':
                    self.rxn_list[j].initial = []
                    LatState = RawTxt[i+1:i+1+nSites]
                    for k in range(0, len(LatState)):
                            self.rxn_list[j].initial.\
                            append(LatState[k].split('\n')[0])
                    for k in range(0, nSites):
                        StateLine.append(i+1+k)
                elif RawTxt[i].split()[0] == 'final':
                    self.rxn_list[j].final = []
                    LatState = RawTxt[i + 1:i + 1 + nSites]
                    for k in range(0, len(LatState)):
                            self.rxn_list[j].final.\
                            append(LatState[k].split('\n')[0])
                    for k in range(0, nSites):
                        StateLine.append(i+1+k)

            if n_variants == 0:     # There are no variants, just one version of the reaction
            
                #raise NameError('No variants for reaction: ' + self.rxn_list[j].name)
                n_variants = 1
                
                # Find the beginning and end of the reaction information
                variantInd = np.array([[0, 0]]*n_variants)
                variantInd[0, 0] = MechInd[j, 0] + 1 
                variantInd[0, 1] = MechInd[j, 1]
                
            else:
                    
                # Find beginning and ending lines for each variant
                variantInd = np.array([[0, 0]]*n_variants)
                Count  = 0
                for i in range(MechInd[j, 0] + 1, MechInd[j, 1]):
                    if RawTxt[i].split()[0] == 'variant':
                        variantInd[Count, 0] = i
                    if RawTxt[i].split()[0] == 'end_variant':
                        variantInd[Count, 1] = i
                        Count += 1

            self.rxn_list[j].variant_list = [ rxn_variant() for i in range(n_variants) ]
                    
            # Loop over list of recation variants        
            for k in range( n_variants ):
            
                for i in range(variantInd[k, 0], variantInd[k, 1]):
                    if RawTxt[i].split()[0] == 'variant':
                        self.rxn_list[j].variant_list[k].name = RawTxt[i].split()[1]
                    elif RawTxt[i].split()[0] == 'site_types':
                        self.rxn_list[j].variant_list[k].site_types = RawTxt[i].split()[1:]
                    elif RawTxt[i].split()[0] == 'pre_expon':
                        self.rxn_list[j].variant_list[k].pre_expon = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'pe_ratio':
                        self.rxn_list[j].variant_list[k].pe_ratio = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'activ_eng':
                        self.rxn_list[j].variant_list[k].activ_eng = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'prox_factor':
                        self.rxn_list[j].variant_list[k].prox_factor = float(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == '#':
                        pass
                
                # If there is no variant, it does not have an extra name
                if self.rxn_list[j].variant_list[k].name is None:
                    self.rxn_list[j].variant_list[k].name = ''
                
                # Assign scaledown factor if it is present
                if StiffCorrLine != -1:
                    self.rxn_list[j].variant_list[k].scaledown_factor = scaledown_factor_list[all_rxn_ind]
                all_rxn_ind += 1


    def WriteIn(self, fldr):
    
        '''
        Write mechanism_input.dat
        '''

        with open(os.path.join(fldr, self.fname), 'w') as txt:
        
            txt.write('mechanism\n\n')
            
            if self.include_scaledown:
                txt.write('# Automated stiffness reconditioning employed\n')
                txt.write('# \n')
                txt.write('# SDF: ')
                for i in self.rxn_list:
                    for j in i.variant_list:
                        txt.write('{0:.5e} \t'.format( j.scaledown_factor ))
                txt.write('\n\n')
            
            # Loop through reactions            
            for rxn in self.rxn_list :

                txt.write('#'*80 + '\n\n')

                if rxn.is_reversible:
                    txt.write('reversible_step ' + rxn.name + '\n')
                else:
                    txt.write('step ' + rxn.name + '\n')

                txt.write('  sites ' + str(rxn.sites) + '\n')
                if not rxn.neighboring is None:
                    txt.write('  neighboring')
                    for j in range(0, len(rxn.neighboring)):
                        txt.write(' ' + rxn.neighboring[j])
                    txt.write('\n')

                if not rxn.gas_reacs_prods is None:
                    txt.write('  gas_reacs_prods')
                    for j in range(0,len(rxn.gas_reacs_prods)):
						txt.write(' ' + str(rxn.gas_reacs_prods[j]))
                    txt.write('\n')

                txt.write('  initial\n')
				
                for j in range( rxn.sites ):
                    txt.write(rxn.initial[j] + '\n')

                txt.write('  final\n')
                for j in range( rxn.sites ):
                    txt.write(rxn.final[j] + '\n')

                txt.write('\n')
                
                # Loop through variants
                for rxn_var in  rxn.variant_list :
                    txt.write('  {} {}\n'.format('variant', rxn_var.name))
                    txt.write('    {:25}'.format('site_types'))
                    for k in range( len( rxn_var.site_types )):
                        txt.write('{} '.format( rxn_var.site_types[k]) )
                    txt.write('\n')
                    
                    # Write pre-exponential factor. Add comment if it has been rescaled
                    if rxn_var.scaledown_factor == 1.0:             # reaction has not been rescaled
                        txt.write('    {:25}{:.5e}\n'.format('pre_expon', rxn_var.pre_expon))
                    else:                                   # reaction has been rescaled
                        txt.write('    {:25}{:.5e}'.format('pre_expon', rxn_var.pre_expon))
                        txt.write( ('    # Pre-exponential has been ' + 'rescaled by a factor of {0:.5e}\n').format(rxn_var.scaledown_factor) )

                    if rxn.is_reversible:
                        txt.write('    {:25}{:.5e}\n'.format('pe_ratio', rxn_var.pe_ratio))
                        txt.write('    {:25}{:5.3f}\n'.format('prox_factor', (rxn_var.prox_factor) ) )
                        
                    txt.write('    {:25}{:4.2f}\n'.format('activ_eng', (rxn_var.activ_eng )) )
                        
                    txt.write('  end_variant\n\n')

                if rxn.is_reversible:
                    txt.write('end_reversible_step\n\n')
                else:
                    txt.write('end_step\n\n')

            txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_mechanism')

    def CalcThermo(self, fldr, T):
        '''
        Calculate the forward activation energy, forward and reverse
        pre-exponential factors and the PE-ratio for each reaction described
        in Mechanism_input.dat using an input file with energies and
        vibrational frequencies for all species and transition states

        Assumes that each reaction has only 1 variant (MPN)
        '''
        filepath = os.path.join(fldr, 'Zacros_Species_Energy.txt')
        species_data = read_excel(io=filepath)
        T_species = [Zacros(**specie_data) for specie_data in species_data]

        '''
        Create list of transition state species
        '''
        TST = []

        for y in range(0, len(T_species)):
            if T_species[y].name.startswith('TST') or T_species[y].name == '*':
                TST.append([T_species[y].name, y])

        '''
        Recalculate all entries in mechanism_input.dat
        '''
        for x in range(0, len(self.rxn_list)):
            q_vib_surf = []
            Rxn_TST = 'TST' + ('0' + str(x + 1))[-2:]
            TST_index = -1
            '''
            Find the index of the transition state species and the slab
            energy species for the current reaction
            '''
            for element in TST:
                if element[0] == Rxn_TST:
                    TST_index = element[1]
                elif element[0] == '*':
                    TST_Slab = element[1]
            '''
            Create list of all surface products and reactants
            '''
            surf_species = []
            for e in self.rxn_list[x].initial:
                surf_species.append(e.split())
            surf_prod = []
            for e in self.rxn_list[x].final:
                surf_prod.append(e.split())
            activ_eng = 0.0
            fwd_pre = 0.0

            if TST_index == -1:
                '''
                Case = No transition state energetics provided
                '''
                if not self.rxn_list[x].gas_reacs_prods is None:
                    MW_gas = next(e.MW for e in T_species
                                  if e.name == self.rxn_list[x].
                                  gas_reacs_prods[0])
                    q_vib_gas = next(e.q_vib for e in T_species
                                     if e.name == self.rxn_list[x].
                                     gas_reacs_prods[0])
                    q_rot_gas = next(e.q_rot for e in T_species
                                     if e.name == self.rxn_list[x].
                                     gas_reacs_prods[0])
                    q_trans2D_gas = next(e.q_trans2D for e in T_species
                                         if e.name == self.rxn_list[x].
                                         gas_reacs_prods[0])
                    for y in range(0, len(surf_prod)):
                        if surf_prod[y][1] != '*' and\
                                              int(surf_prod[y][2]) == 1:
                            q_vib_surf.append(next(e.q_vib for e in T_species
                                                   if e.name ==
                                                   surf_prod[y][1]))

                if not self.rxn_list[x].gas_reacs_prods is None and\
                   int(self.rxn_list[x].gas_reacs_prods[1]) == -1:
                    '''
                    No transition state and a gas reactant
                    Non-activated adsorbtion
                    '''
                    fwd_pre = T_species[x].A_st /\
                        np.sqrt(2*np.pi * MW_gas * _c.kb1*T)\
                        * 1e5

                    rev_pre = q_vib_gas * q_rot_gas * q_trans2D_gas /\
                        np.product(q_vib_surf) * _c.kb1 * T/_c.h1

                elif not self.rxn_list[x].gas_reacs_prods is None and\
                        int(self.rxn_list[x].gas_reacs_prods[1]) == 1:
                    '''
                    No transition state and a gas product
                    Non-activated desorbtion
                    '''
                    rev_pre = T_species[x].A_st /\
                        np.sqrt(2*np.pi * MW_gas * _c.kb1*T)\
                        * 1e5

                    fwd_pre = q_vib_gas * q_rot_gas * q_trans2D_gas /\
                        np.product(q_vib_surf) *\
                        _c.kb1 * T/_c.h1
                else:
                    '''
                    Insufficient information to calculate pre-exponential
                    factors.  Set values to zero
                    '''
                    fwd_pre = 0
                    rev_pre = 1
            else:
                '''
                Case = Transition state energetics provided
                '''
                activ_eng = T_species[TST_index].etotal -\
                    T_species[TST_Slab].etotal +\
                    T_species[TST_index].zpe/_c.ev_atom_2_kcal_mol
                q_vib_TST = next(e.q_vib for e in T_species
                                 if e.name ==
                                 T_species[TST_index].name)
                if not self.rxn_list[x].gas_reacs_prods is None:
                    q_vib_gas = next(e.q_vib for e in T_species
                                     if e.name == self.rxn_list[x].
                                     gas_reacs_prods[0])
                    q_rot_gas = next(e.q_rot for e in T_species
                                     if e.name == self.rxn_list[x].
                                     gas_reacs_prods[0])
                    q_trans2D_gas = next(e.q_trans2D for e in T_species
                                         if e.name == self.rxn_list[x].
                                         gas_reacs_prods[0])
                    A_st = next(e.A_st for e in T_species
                                if e.name ==
                                self.rxn_list[x].gas_reacs_prods[0])
                    MW_gas = next(e.MW for e in T_species
                                  if e.name ==
                                  self.rxn_list[x].gas_reacs_prods[0])
                    for y in range(0, len(surf_prod)):
                        if surf_prod[y][1] != '*' and\
                              int(surf_prod[y][2]) == 1:
                                q_vib_surf.append(next(e.q_vib
                                                  for e in T_species
                                                  if e.name ==
                                                  surf_prod[y][1]))
                    Q_gas = q_vib_gas * q_rot_gas * q_trans2D_gas

                if not self.rxn_list[x].gas_reacs_prods is None and\
                   int(self.rxn_list[x].gas_reacs_prods[1]) == -1:
                    '''
                    Transition state and a gas reactant
                    Activated adsorbtion
                    '''
                    activ_eng -=\
                        next(e.etotal + e.zpe/_c.ev_atom_2_kcal_mol
                             for e in T_species
                             if e.name == self.rxn_list[x].gas_reacs_prods[0])
                    fwd_pre = q_vib_TST/Q_gas * A_st /\
                        np.sqrt(2*np.pi*MW_gas*_c.kb1*T)*1e5
                    rev_pre = q_vib_TST/np.product(q_vib_surf) *\
                        (_c.kb1*T/_c.h1)
                elif not self.rxn_list[x].gas_reacs_prods is None and\
                        int(self.rxn_list[x].gas_reacs_prods[1]) == 1:
                    '''
                    Transition state and a gas product
                    Activated desorbtion
                    '''
                    q_vib_surf = []
                    for y in range(0, len(surf_species)):
                        if surf_species[y][1] != '*' and\
                              int(surf_species[y][2]) == 1:
                                activ_eng -=\
                                    (next(e.etotal + e.zpe /
                                          _c.ev_atom_2_kcal_mol
                                          for e in T_species
                                     if e.name == surf_species[y][1]) -
                                     T_species[TST_Slab].etotal)
                                q_vib_surf.append(next(e.q_vib
                                                  for e in T_species
                                                  if e.name ==
                                                  surf_species[y][1]))
                    rev_pre = q_vib_TST/Q_gas * A_st /\
                        np.sqrt(2*np.pi*MW_gas*_c.kb1*T)*1e5
                    fwd_pre = q_vib_TST/np.product(q_vib_surf) *\
                        (_c.kb1*T/_c.h1)
                else:
                    '''
                    Transition state and no gas reactant or product
                    Surface reaction
                    '''
                    q_vib_surf = []
                    for y in range(0, len(surf_species)):
                        if int(surf_species[y][2]) == 1 and\
                          surf_species[y][1] != '*':
                            activ_eng -=\
                                (next(e.etotal + e.zpe/_c.ev_atom_2_kcal_mol
                                      for e in T_species
                                 if e.name == surf_species[y][1]) -
                                 T_species[TST_Slab].etotal)
                            q_vib_surf.append(next(e.q_vib for e in T_species
                                              if e.name == surf_species[y][1]))
                    q_vib_reactants = np.product(q_vib_surf)
                    fwd_pre = q_vib_TST/q_vib_reactants * (_c.kb1*T/_c.h1)

                    q_vib_prod = []
                    for y in range(0, len(surf_prod)):
                        if int(surf_prod[y][2]) == 1 and\
                          surf_prod[y][1] != '*':
                            q_vib_prod.append(next(e.q_vib for e in T_species
                                                   if e.name ==
                                                   surf_prod[y][1]))
                    q_vib_products = np.product(q_vib_prod)
                    rev_pre = q_vib_TST/q_vib_products * (_c.kb1*T/_c.h1)
                    
            # Modify reaction data
            # assumes that each reaction has only one variant - i.e. variant_list[0]
            self.rxn_list[x].variant_list[0].activ_eng = max(activ_eng, 0.0)
            self.rxn_list[x].variant_list[0].pre_expon = fwd_pre *\
                self.rxn_list[x].variant_list[0].scaledown_factor
            self.rxn_list[x].variant_list[0].pe_ratio = fwd_pre/rev_pre
        




class StateIn(object):

    '''
    Handles input from state_input.dat
    '''
    
    fname = 'state_input.dat'

    def __init__(self):
    
        self.Type = None    # None: No state_input.dat, StateInput: has read state_input.dat, history: read from history file previously
        self.Struct = None
        
    def ReadIn(self, fldr):
    
        '''
        Read state_input.dat
        '''

        if os.path.isfile(os.path.join(fldr, 'state_input.dat')):
        
            with open(os.path.join(fldr, 'state_input.dat'), 'r') as Txt:
                RawTxt = Txt.readlines()
                
            self.Struct = []
            for i in RawTxt:
                self.Struct.append(i.split('\n')[0])
                
            self.Type = 'StateInput'
            
        else:
            self.Type = None
        

    def WriteIn(self, fldr, surf_spec):
    
        '''
        Write state_input.dat
        '''

        if self.Type is None:
            pass

        elif self.Type == 'StateInput':
            with open(os.path.join(fldr, self.fname), 'w') as txt:
                for i in self.Struct:
                    txt.write(i + '\n')

        elif self.Type == 'history':

            Lattice = self.Struct
            UniqSpec = np.unique(Lattice[np.not_equal(
                    Lattice[:, 2], 0), 1])
            nAds = len(UniqSpec)
            SpecIden = [0] * nAds
            AdsInfo = [[] for i in range(0, nAds)]
            DentInfo = [[] for i in range(0, nAds)]
            for i in range(0, nAds):
                for j in range(0, Lattice.shape[0]):
                    if UniqSpec[i] == Lattice[j, 1]:
                        AdsInfo[i].append(j + 1)
                        DentInfo[i].append(Lattice[j, 3])
                        SpecIden[i] = Lattice[j, 2]

            if nAds > 0:
                with open(os.path.join(fldr, self.fname), 'w') as txt:
                    txt.write('initial_state\n')
                    for i in range(0, nAds):
                        txt.write('  seed_on_sites  {:10}'.
                                  format(surf_spec
                                         [SpecIden[i]-1], 10))
                        for j in range(0, len(DentInfo[i])):
                            for k in range(0, len(DentInfo[i])):
                                if j + 1 == DentInfo[i][k]:
                                    txt.write(str(AdsInfo[i][k]) + '  ')
                        txt.write('\n')
                    txt.write('end_initial_state\n')
        else:
            print 'Unrecognized state_input type'
            print 'state_input not written'


'''
============ Classes to handle output files ============
'''
            
class PerformanceOut(object):

    '''
    Handles data from general_output.txt
    '''
    
    fname = 'general_output.txt'

    def __init__(self):
        
        self.nRxn = None
        self.t_final = None
        self.events_occurred = None
        self.CPU_time = None
        self.RxnNameList = []
        self.Nu = None
        self.UniqNu = None
        
    
    def ReadOut(self, fldr, surf_spec_names, gas_spec_names):
        '''
        Read general_output.txt
        '''
        
        n_surf = len(surf_spec_names)
        n_gas = len(gas_spec_names)
        
        with open(os.path.join(fldr, self.fname), 'r') as txt:
            RawTxt = txt.readlines()
    
        for i in range(0, len(RawTxt)):
            if _re.search('Number of elementary steps:', RawTxt[i]):
                self.nRxn = np.int(RawTxt[i].split(':')[1])
            elif _re.search('Current KMC time:', RawTxt[i]):
                self.t_final = np.float(RawTxt[i].split(':')[1])
            elif _re.search('Events occurred:', RawTxt[i]):
                self.events_occurred =\
                np.float(RawTxt[i].split(':')[1])
            elif _re.search('Elapsed CPU time:', RawTxt[i]):
                after_colon = RawTxt[i].split(':')[1]
                self.CPU_time =\
                    np.float(after_colon.split(' ')[-2])
            elif _re.search('Reaction network:', RawTxt[i]):
                RxnStartLine = i + 2
    
        if RawTxt[RxnStartLine].split()[0] == '1.':
            NameInd = 1
        else:
            NameInd = 0
    
        self.RxnNameList = []
        nuList = []
        for i in range(RxnStartLine, RxnStartLine + self.nRxn):
            RxnName = RawTxt[i].split()[NameInd][:-1]
            self.RxnNameList.append(RxnName)
            RxnStr = RawTxt[i][_re.search('Reaction:', RawTxt[i]).end():]
            RxnStrList = RxnStr.split()
            nu = [0] * (n_surf + n_gas)
            for j in range(0, len(RxnStrList)):
                if RxnStrList[j] == '->':
                    ArrowInd = j
            for j in range(0, len(RxnStrList)):
                if j < ArrowInd:
                    Sign = -1
                else:
                    Sign = 1
    
                if _re.search('\(', RxnStrList[j]):
                    SurfIden = _re.sub(r'\([^)]*\)', '', RxnStrList[j])
                    if SurfIden != '*':
                        SurfInd = [k for k in
                                range( n_surf )
                                if SurfIden ==
                                surf_spec_names[k]][0]
                        nu[SurfInd] += Sign
                elif RxnStrList[j] != '->' and RxnStrList[j] != '+':
                    GasInd = [k for k in
                            range( n_gas )
                            if RxnStrList[j] ==
                            gas_spec_names[k]][0]
                    nu[ n_surf + GasInd] += Sign
            nuList.append(nu)
    
        self.Nu = nuList
        self.UniqNu = _ReturnUnique(nuList).tolist()


class ProcstatOut(object):
    
    '''
    Handles data from procstat_output.txt
    '''
    
    fname = 'procstat_output.txt'

    def __init__(self):
    
        self.Spacing = None
        self.t = None
        self.events = None
    
    def ReadOut(self, fldr):
        '''
        Read procstat_output.txt
        '''
        MaxLen = np.int(2e4)
        with open(os.path.join(fldr, self.fname), 'r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen * 3:  # Procstat uses 3 lines per outputs
            Spacing = np.int(np.floor((len(RawTxt) - 1)/(MaxLen*3)))
            RawTxt2 = []
            for i in range(0, MaxLen):
                RawTxt2.append(RawTxt[i*Spacing*3+1])
                RawTxt2.append(RawTxt[i*Spacing*3+2])
                RawTxt2.append(RawTxt[i*Spacing*3+3])
        else:
            Spacing = 1
            RawTxt2 = RawTxt[1:]

        t = []
        events = []
        for i in range(0, len(RawTxt2)/3):
            t.append(np.float(RawTxt2[i*3].split()[3]))
            eventsTemp = RawTxt2[i*3+2].split()[1:]
            for j in range(0, len(eventsTemp)):
                eventsTemp[j] = np.int(eventsTemp[j])
            events.append(eventsTemp)

        self.Spacing = Spacing
        self.t = np.asarray(t)
        self.events = np.asarray(events)


class SpecnumOut(object):
    
    '''
    Handles data from specnum_output.txt
    '''
    
    fname = 'specnum_output.txt'
    
    def __init__(self):
        
        '''
        Initializes class variables
        '''
    
        self.Spacing = None
        self.nEvents = None
        self.t = None
        self.T = None
        self.E = None
        self.spec = None
        
    
    def ReadOut(self, fldr):
        '''
        Read specnum_output.txt
        '''
        MaxLen = np.int(2e4)
        with open(os.path.join(fldr, self.fname), 'r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen:
            Spacing = np.int(np.floor((len(RawTxt)-1)/MaxLen))
            RawTxt2 = []
            for i in range(0, MaxLen):
                RawTxt2.append(RawTxt[i*Spacing+1])
        else:
            Spacing = 1
            RawTxt2 = RawTxt[1:]

        nEvents = []
        t = []
        T = []
        E = []
        spec = []

        for i in range(0, len(RawTxt2)):
            LineSplit = RawTxt2[i].split()
            nEvents.append(np.int(LineSplit[1]))
            t.append(np.float(LineSplit[2]))
            T.append(np.float(LineSplit[3]))
            E.append(np.float(LineSplit[4]))
            specTemp = LineSplit[5:]
            for j in range(0, len(specTemp)):
                specTemp[j] = np.int(specTemp[j])
            spec.append(specTemp)
            
        # Store data in class variables
        self.Spacing = Spacing
        self.nEvents = np.asarray(nEvents)
        self.t = np.asarray(t)
        self.T = np.asarray(T)
        self.E = np.asarray(E)
        self.spec = np.asarray(spec)
    

class HistoryOut():

    '''
    Handles data from history_output.txt
    '''

    fname = 'history_output.txt'
    
    def __init__(self):
    
        n_snapshots = 0
        snapshots = []
        snap_times = None
        
    def ReadOut(self, fldr, nSites):
    
        '''
        Read history_output.txt
        fldr: name of the folder containting the file
        nSites: number of lattice sites, obtained from lattice_output.txt
        '''
        
        HistPath = os.path.join(fldr, self.fname)
        
        # Check if file exists
        if not os.path.isfile(HistPath):
            return

        nLines = _rawbigcount(HistPath)
        
        self.n_snapshots = (nLines-6)/(nSites+2)
        self.snapshots = []
        self.snap_times = []

        for snap_ind in range(self.n_snapshots):
            snap_data = np.array([[0]*4]*nSites)
            _linecache.clearcache()
            snap_header = _linecache.getline(HistPath, 8 + snap_ind *
                                             (nSites+2)-1).split()
            self.snap_times.append(np.float(snap_header[3]))
            for i in range(0, nSites):
                snap_data[i, :] = _linecache.getline(HistPath, 8 + snap_ind *
                                                     (nSites+2)+i).split()
            self.snapshots.append(snap_data)
        

def StateInc(i):
    if _re.search('off', i):
        state = 'off'
        inc = ''
    elif _re.search('on time', i):
        state = 'time'
        inc = np.float(i.split()[3])
    elif _re.search('on event', i):
        state = 'event'
        if len(i.split()) < 4:
            inc = 1
        else:
            inc = np.int(i.split()[3])
    return (state, inc)
    
    
def Read_propensities(path, nRxn):
    '''
    Read propenisty data from output files. The initial time point: At t = 0 or after 0 events is wrong because
    it is recorded before the first propensities are calculated. Instead, they are erroneously all zeros and should
    not be used for averaging. The rest are correct.        
    
    :param Mode: 0 - Read Prop_output.bin, instantaneous propensities
        1 - Read PropCounter_output.bin, time integrated propensities used for accurate time averages
        
    :returns: Matrix of time integrated surface species populations
    '''
    
    if os.path.isfile(os.path.join(path, 'Prop_output.bin')):
    
        dt = np.dtype(np.float64)
        virtual_arr = np.memmap(os.path.join(path, 'Prop_output.bin'), dt, "r")
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]
    
        prop = np.reshape(virtual_arr, [nNum/nRxn, nRxn])

        del virtual_arr
        return np.array(prop)
        
    elif os.path.isfile(os.path.join(path, 'propensities_output.txt')):
    
        with open(os.path.join(path, 'propensities_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()
        prop = []
        for i in range(len(RawTxt)):
            if i >= 1:
				LineSplit = RawTxt[i].split()
				line_data = []
				if len(LineSplit) > 0:
					for dub in LineSplit:
						line_data.append(np.float(dub))
					prop.append(line_data)
        return np.array(prop)
    
    else:
        return None
    

def Read_time_integrated_propensities(path, nRxn):
    '''
    Read propenisty data from output files. The initial time point: At t = 0 or after 0 events is wrong because
    it is recorded before the first propensities are calculated. Instead, they are erroneously all zeros and should
    not be used for averaging. The rest are correct.        
    
    :param Mode: 0 - Read Prop_output.bin, instantaneous propensities
        1 - Read PropCounter_output.bin, time integrated propensities used for accurate time averages
        
    :returns: Matrix of time integrated surface species populations
    '''
    
    if os.path.isfile(os.path.join(path, 'PropCounter_output.bin')):
    
        dt = np.dtype(np.float64)
        virtual_arr = np.memmap(os.path.join(path, 'PropCounter_output.bin'), dt, "r")
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]
    
        propCounter = np.reshape(virtual_arr, [nNum/nRxn, nRxn])
            
        del virtual_arr
        return np.array(propCounter)
    
    elif os.path.isfile(os.path.join(path, 'timeintprop_output.txt')):
    
        with open(os.path.join(path, 'timeintprop_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()
        propCounter = []
        for i in range(len(RawTxt)):
            
            LineSplit = RawTxt[i].split()
            line_data = []
            if len(LineSplit) > 0:
                for dub in LineSplit:
                    line_data.append(np.float(dub))
                propCounter.append(line_data)
        
        return np.array(propCounter)
    else:    
        return None

    
def Read_trajectory_derivatives(path, nRxn):
    '''
    Read SA_output.bin - get trajectory derivatives for use in likelihood ratio sensitivity analysis
    
    :returns: Matrix of trajectory derivatives
    '''
    dt = np.dtype(np.float64)
    FileName = 'SA_output.bin'
    if os.path.isfile(os.path.join(path, FileName)):
        virtual_arr = np.memmap(os.path.join(path,
                                               FileName), dt, "r")
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]

        W_sen_anal = np.reshape(virtual_arr, [nNum/nRxn, nRxn])

        del virtual_arr
        return np.array(W_sen_anal)
    
    elif os.path.isfile(os.path.join(path, 'trajderiv_output.txt')):
        with open(os.path.join(path, 'trajderiv_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()
        W_sen_anal = []
        for i in range(len(RawTxt)):
            
            LineSplit = RawTxt[i].split()
            line_data = []
            if len(LineSplit) > 0:
                for dub in LineSplit:
                    line_data.append(np.float(dub))
                W_sen_anal.append(line_data)
        return np.array(W_sen_anal)
    else:    
        return None


def Read_time_integrated_species(path, n_surf_specs):
    '''
    Read time integrated species counts
    
    :returns: Matrix of time integrated surface species populations
    '''
    
    if os.path.isfile(os.path.join(path, 'IntegSpec_output.bin')):
        dt = np.dtype(np.float64)
        virtual_arr = np.memmap(os.path.join(path, 'IntegSpec_output.bin'), dt, "r")
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % n_surf_specs)
        virtual_arr = virtual_arr[:nNum]

        spec_num_int = np.reshape(virtual_arr, [nNum/n_surf_specs, n_surf_specs])

        del virtual_arr
        return np.array(spec_num_int)
        
    elif os.path.isfile(os.path.join(path, 'timeintspecs_output.txt')): 
        with open(os.path.join(path, 'timeintspecs_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()
        spec_num_int = []
        for i in range(len(RawTxt)):
            
            LineSplit = RawTxt[i].split()
            line_data = []
            if len(LineSplit) > 0:
                for dub in LineSplit:
                    line_data.append(np.float(dub))
                spec_num_int.append(line_data)
        return np.array(spec_num_int)
        
    else:    
        return None
        
        
def Read_time_integrated_site_props(path, nSites, nRxn, nSnaps ):
    '''
    Read time integrated site propensities
    
    :returns: List of [nSites x nRxn] matrices, one for each snapshot
    '''
    
    
    fname = os.path.join(path, 'TIsiteprops_output.txt')
    if os.path.isfile(fname):
        with open(fname, 'r') as txt:
            RawTxt = txt.readlines()
    
        TS_site_props_list = []
        line_count = 0
        for snap_ind in range(nSnaps):
            
            prop_array = np.zeros([nSites, nRxn])
            for site_ind in range(nSites):
                x = RawTxt[line_count].split()
                for rxn_ind in range(nRxn):
                    prop_array[site_ind, rxn_ind] = x[rxn_ind]
                line_count += 1
                
            line_count += 1
                
            TS_site_props_list.append(prop_array)    
        
        return TS_site_props_list
        
    else:
        return None
            