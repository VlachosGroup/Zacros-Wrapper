# -*- coding: utf-8 -*-
'''
         -----------------------------------------------------
               Read all Zacros input and output files.
            Calculate cluster energies, pre-exponential
             factors and activation energeies using the
           DFT_to_Thermochemistry.py code.  Replace those
           values in the Zacros energetics_input.dat and
            mechanism_input.dat files and write the new
                         file versions.

                     Vlachos Research Group
                Chemical and Biomolecular Egineering
                      University of Delaware

                     Gerhard R Wittreich, P.E.
          -----------------------------------------------------

Created on Fri Apr 28 10:49:53 2017

@author: wittregr

Adopted from Matlab code written and modified by:

                        Marcel Nunez
                            and
                        Taylor Robie

 This program contains the class objects used to read energy, vibration and
 molecular configuration data and determine the standard entropy and enthalpy
 and heat capacities at various temperatures.

'''

import os as _os
import numpy as _np
import re as _re
import random as _random
import linecache as _linecache
import DFT_to_Thermochemistry as _thermo
import IOdata_v2 as _ClassIn

from GRW_constants import constant as _c
from Helper import ReadWithoutBlankLines as _ReadWithoutBlankLines
from Helper import ReturnUnique as _ReturnUnique
from Helper import rawbigcount as _rawbigcount
reload(_thermo)


class IOdata(object):

    def __init__(self):
        self.Path = 'C:\Users\wittregr\Documents\Python Scripts'
        self.Cluster_Num = 0
        self.Reaction_Num = 0
        self.C_index = 0
        self.R_index = 0
        self.V_index = 0

    def ReadAllInput(self):
        '''
        Read all input files. state_input.dat will be read only
        if it is present.
        '''
        self.ReadSimIn()
        self.ReadLatticeIn()
        self.ReadEngIn()
        self.ReadMechIn()

        if _os.path.isfile(_os.path.join(self.Path, 'state_input.dat')):
            self.ReadStateInput()
        else:
            self.State = None

    def WriteAllInput(self):
        '''
        Write all Zacros input files
        '''
        self.WriteSimIn()
        self.WriteMechanism()
        self.WriteEnergetics()
        self.WriteStateIn()
        self.WriteLattice()

    def ReadAllOutput(self, build_lattice=False):
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
                self.KMC_lat.Read_lattice_output(_os.path.join
                                                 (self.Path,
                                                  'lattice_output.txt'))
            self.ReadHistory()

            # Extra binary files
            if _os.path.isfile(_os.path.join(self.Path,
                                             'PropCounter_output.bin')):
                self.ReadProp(1)
            if _os.path.isfile(_os.path.join(self.Path,
                                             'SA_output.bin')):
                self.ReadSA()
#            self.ReadCluster()
#            self.ReadProp(0)

        else:
            print 'general_output.txt not found in ' + self.Path

    def FindCluster(self, Cluster_Num):
        '''
        Method finds the Cluster and Variant index of the nth
        Cluster-Variant where n is specified by Cluster_Num and
        the indices are returned as C.index (Cluster)
        and V.index (Variant) such that Cluster[C_index].variant_name[V_index]
        represents the name of the nth Cluster-Variant
        '''
        Cluster_Num = int(Cluster_Num)
        Tvariants = sum(s.nVariant for s in self.Cluster)
        if Tvariants >= Cluster_Num and Cluster_Num >= 1:
            var = []
            for s in self.Cluster:
                var.append(s.nVariant)
            var = _np.array(var)
            C_index = _np.argmin(var.cumsum() < Cluster_Num)
            V_index = Cluster_Num - sum(var[0:C_index]) - 1
        else:
            C_index = V_index = -1
        return(C_index, V_index)

    def FindReaction(self, Reaction_Num):
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
            var = _np.array(var)
            R_index = _np.argmin(var.cumsum() < Reaction_Num)
            V_index = Reaction_Num - sum(var[0:R_index]) - 1
        else:
            R_index = V_index = -1
        return(R_index, V_index)

    def CheckComplete(self):
        '''
        Check to see if a Zacros run has completed successfully
        '''
        Complete = False
        if _os.path.isfile(_os.path.join(self.Path,
                                         'general_output.txt')):
            with open(_os.path.join(self.Path,
                                    'general_output.txt'),
                      'r') as txt:
                RawTxt = txt.readlines()
            for i in RawTxt:
                if _re.search('Normal termination', i):
                    Complete = True
        return Complete

    def ReadEngIn(self):
        '''
        Read energetics_input.dat
        '''
        RawTxt = _ReadWithoutBlankLines(_os.path.join(self.Path,
                                                      'energetics_input.dat'),
                                        CommentLines=False)
        nLines = len(RawTxt)

        nCluster = 0
        for i in range(0, nLines):
            if RawTxt[i].split()[0] == 'cluster':
                nCluster += 1

        ClusterInd = _np.array([[0, 0]]*nCluster)
        Count = 0
        for i in range(0, nLines):
            if RawTxt[i].split()[0] == 'cluster':
                ClusterInd[Count, 0] = i
            if RawTxt[i].split()[0] == 'end_cluster':
                ClusterInd[Count, 1] = i
                Count += 1

        nClusterTotal = 0
        self.Cluster = []
        for j in range(0, nCluster):
            self.Cluster.append(_ClassIn.ClusterIn())
            self.Cluster[j].name = RawTxt[ClusterInd[j, 0]].split()[1]
            Count = 0
            for i in range(ClusterInd[j, 0] + 1, ClusterInd[j, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    Count += 1
                elif RawTxt[i].split()[0] == 'sites':
                    self.Cluster[j].sites = int(RawTxt[i].split()[1])
                elif RawTxt[i].split()[0] == 'neighboring':
                    self.Cluster[j].neighboring = RawTxt[i].split()[1:]
                elif RawTxt[i].split()[0] == 'lattice_state':
                    self.Cluster[j].latstate = RawTxt[i + 1:i + 1 +
                                                      self.Cluster[j].sites]
                    for k in range(0, len(self.Cluster[j].latstate)):
                        self.Cluster[j].latstate[k] =\
                         self.Cluster[j].latstate[k].split('\n')[0]

            self.Cluster[j].nVariant = Count
            nClusterTotal += self.Cluster[j].nVariant
            variantInd = _np.array([[0, 0]]*self.Cluster[j].nVariant)
            Count = 0
            for i in range(ClusterInd[j, 0]+1, ClusterInd[j, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    variantInd[Count, 0] = i
                if RawTxt[i].split()[0] == 'end_variant':
                    variantInd[Count, 1] = i
                    Count += 1

            for k in range(0, self.Cluster[j].nVariant):
                for i in range(variantInd[k, 0], variantInd[k, 1]):
                    if RawTxt[i].split()[0] == 'variant':
                        self.Cluster[j].variant_name.\
                         append(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'site_types':
                        self.Cluster[j].variant_site_types.\
                         append(RawTxt[i].split()[1:])
                    elif RawTxt[i].split()[0] == 'graph_multiplicity':
                        self.Cluster[j].variant_graph_multiplicity.\
                         append(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'cluster_eng':
                        self.Cluster[j].variant_eng.\
                         append(float(RawTxt[i].split()[1]))

    def WriteEnergetics(self):
        '''
        Write energetics_input.dat
        '''
        if not hasattr(self, 'Cluster'):
            self.ReadEngIn()

        nCluster = len(self.Cluster)

        with open(_os.path.join(self.Path,
                                'energetics_input.dat'), 'w') as txt:
            txt.write('energetics\n\n')
            for i in range(0, nCluster):
                txt.write('#'*80 + '\n\n')
                txt.write('cluster ' + self.Cluster[i].name + '\n\n')
                txt.write('  sites ' + str(self.Cluster[i].sites) + '\n')

                if hasattr(self.Cluster[i], 'neighboring'):
                    txt.write('  neighboring')
                    for j in range(0, len(self.Cluster[i].neighboring)):
                        txt.write(' ' + self.Cluster[i].neighboring[j])
                    txt.write('\n')

                txt.write('  lattice_state\n')
                for j in range(0, int(self.Cluster[i].sites)):
                    txt.write(self.Cluster[i].latstate[j] + '\n')

                nVariant = len(self.Cluster[i].variant_name)
                txt.write('\n')
                for j in range(0, nVariant):
                    txt.write('  {} {}\n'.format('variant',
                              self.Cluster[i].variant_name[j]))
                    txt.write('    {:25}'.format('site_types'))
                    for k in range(0,
                                   len(self.Cluster[i].variant_site_types[j])):
                        txt.write('{} '.format
                                  (self.Cluster[i].variant_site_types[j][k]))
                    txt.write('\n')
                    if int(self.Cluster[i].variant_graph_multiplicity[j]) > 0:
                        txt.write('    {:25}{}\n'.format('graph_multiplicity',
                                  str(self.Cluster[i].
                                      variant_graph_multiplicity[j])))
                    txt.write('    {:25}{}\n'.format('cluster_eng',
                              str(self.Cluster[i].variant_eng[j])))
                    txt.write('  end_variant\n\n')

                txt.write('end_cluster\n\n')
            txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_energetics')

    def ReadMechIn(self):
        '''
        Read mechanism_input.dat
        '''
        RawTxt = _ReadWithoutBlankLines(_os.path.join(self.Path,
                                                      'mechanism_input.dat'),
                                        CommentLines=True)
        nLines = len(RawTxt)
        StiffCorrLine = -1

        nMech = 0
        for i in range(0, nLines):
            if RawTxt[i].split()[0] == 'reversible_step' or\
                                       RawTxt[i].split()[0] == 'step':
                nMech += 1
            elif _re.search('# Automated stiffness reconditioning employed',
                            RawTxt[i]):
                StiffCorrLine = i

        if StiffCorrLine != -1:
            self.scaledown_factors = [_np.float(i)
                                      for i in RawTxt[StiffCorrLine+2].
                                      split(':')[1].split()]

        # Initialize varibles
        Rxn = []

        MechInd = _np.array([[0, 0]]*nMech)
        Count = 0
        for i in range(0, nLines):

            if RawTxt[i].split()[0] == 'reversible_step':
                MechInd[Count, 0] = i
                Rxn.append(True)

            elif RawTxt[i].split()[0] == 'step':
                MechInd[Count, 0] = i
                Rxn.append(False)

            elif RawTxt[i].split()[0] == 'end_reversible_step':
                MechInd[Count, 1] = i
                Count += 1

            elif RawTxt[i].split()[0] == 'end_step':
                MechInd[Count, 1] = i
                Count += 1

        self.Reaction = []
        for j in range(0, nMech):

            self.Reaction.append(_ClassIn.MechanismIn())
            self.Reaction[j].reversible = Rxn[j]

            # Count the variants

            self.Reaction[j].name = RawTxt[MechInd[j, 0]].split()[1]
            Count = 0
            InVariant = False
            StateLine = []
            for i in range(MechInd[j, 0] + 1, MechInd[j, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    Count += 1
                    InVariant = True
                elif RawTxt[i].split()[0] == 'end_variant':
                    InVariant = False
                elif RawTxt[i].split()[0] == 'gas_reacs_prods':
                    self.Reaction[j].gas_reacs_prods = RawTxt[i].split()[1:]
                elif RawTxt[i].split()[0] == 'sites':
                    nSites = int(RawTxt[i].split()[1])
                    self.Reaction[j].sites = nSites
                elif RawTxt[i].split()[0] == 'neighboring':
                    neighbor = RawTxt[i].split()[1:]
                    self.Reaction[j].neighboring = neighbor
                elif RawTxt[i].split()[0] == 'initial':
                    LatState = RawTxt[i+1:i+1+nSites]
                    for k in range(0, len(LatState)):
                            self.Reaction[j].initial.\
                             append(LatState[k].split('\n')[0])
                    for k in range(0, nSites):
                        StateLine.append(i+1+k)
                elif RawTxt[i].split()[0] == 'final':
                    LatState = RawTxt[i + 1:i + 1 + nSites]
                    for k in range(0, len(LatState)):
                            self.Reaction[j].final.\
                             append(LatState[k].split('\n')[0])
                    for k in range(0, nSites):
                        StateLine.append(i+1+k)
                elif not InVariant and i not in StateLine:
                    print 'Unparsed line in mechanism input:'
                    print RawTxt[i]
            self.Reaction[j].nVariant = Count
            variantInd = _np.array([[0, 0]]*self.Reaction[j].nVariant)
            Count = 0
            for i in range(MechInd[j, 0] + 1, MechInd[j, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    variantInd[Count, 0] = i
                if RawTxt[i].split()[0] == 'end_variant':
                    variantInd[Count, 1] = i
                    Count += 1

            # Enumerate reversibilities

    #        if reversible:
    #            Reactions['nrxns_total'] += 2 * nVariant
    #        else:
    #            Reactions['nrxns_total'] += nVariant

            for k in range(0, self.Reaction[j].nVariant):
                for i in range(variantInd[k, 0], variantInd[k, 1]):
                    if RawTxt[i].split()[0] == 'variant':
                        self.Reaction[j].variant_name.\
                         append(RawTxt[i].split()[1])
                    elif RawTxt[i].split()[0] == 'site_types':
                        self.Reaction[j].variant_site_types.\
                         append(RawTxt[i].split()[1:])
                    elif RawTxt[i].split()[0] == 'pre_expon':
                        self.Reaction[j].variant_pre_expon.\
                         append(float(RawTxt[i].split()[1]))
                    elif RawTxt[i].split()[0] == 'pe_ratio':
                        self.Reaction[j].variant_pe_ratio.\
                         append(float(RawTxt[i].split()[1]))
                    elif RawTxt[i].split()[0] == 'activ_eng':
                        self.Reaction[j].activ_eng.\
                         append(float(RawTxt[i].split()[1]))
                    elif RawTxt[i].split()[0] == 'prox_factor':
                        self.Reaction[j].variant_prox_factor.\
                         append(float(RawTxt[i].split()[1]))
                    elif RawTxt[i].split()[0] == '#':
                        pass
    #                    else:
    #                        print 'Unparsed line in mechanism variant:'
    #                        print RawTxt[i]
        pass
        if StiffCorrLine == -1:
            self.scaledown_factors = _np.ones(sum(len(s.variant_name)
                                              for s in self.Reaction))

    def WriteMechanism(self):
        '''
        Write mechanism_input.dat
        '''
        if not hasattr(self, 'Reaction') or\
           not hasattr(self, 'scaledown_factors'):
            self.ReadMechIn()

        if self.scaledown_factors is None:
            SDBool = False
        else:
            SDBool = True
        nMech = len(self.Reaction)
        StiffCorrCounter = -1
        with open(_os.path.join(self.Path,
                                'mechanism_input.dat'), 'w') as txt:
            txt.write('mechanism\n\n')
            if SDBool:
                txt.write('# Automated stiffness reconditioning employed\n')
                txt.write('# \n')
                txt.write('# SDF: ')
                for i in self.scaledown_factors:
                    txt.write('{0:.5e} \t'.format(i))
                txt.write('\n\n')
            for i in range(0, nMech):

                txt.write('#'*80 + '\n\n')

                if self.Reaction[i].reversible:
                    txt.write('reversible_step ' +
                              self.Reaction[i].name + '\n')
                else:
                    txt.write('step ' + self.Reaction[i].name + '\n')

                txt.write('  sites ' + str(self.Reaction[i].sites) + '\n')
                if hasattr(self.Reaction[i], 'neighboring'):
                    txt.write('  neighboring')
                    for j in range(0, len(self.Reaction[i].neighboring)):
                        txt.write(' ' + self.Reaction[i].neighboring[j])
                    txt.write('\n')

                if hasattr(self.Reaction[i], 'gas_reacs_prods'):
                    txt.write('  {} {} {}\n'.format('gas_reacs_prods',
                              str(self.Reaction[i].gas_reacs_prods[0]),
                              str(self.Reaction[i].gas_reacs_prods[1])))

                txt.write('  initial\n')
                for j in range(0, int(self.Reaction[i].sites)):
                    txt.write(self.Reaction[i].initial[j] + '\n')

                txt.write('  final\n')
                for j in range(0, int(self.Reaction[i].sites)):
                    txt.write(self.Reaction[i].final[j] + '\n')

                nVariant = len(self.Reaction[i].variant_name)
                txt.write('\n')
                for j in range(0, nVariant):
                    txt.write('  {} {}\n'.format('variant',
                              self.Reaction[i].variant_name[j]))
                    txt.write('    {:25}'.format('site_types'))
                    for k in range(0, len(self.Reaction[i].
                                          variant_site_types[j])):
                        txt.write('{} '.format(
                                self.Reaction[i].variant_site_types[j][k]))
                    txt.write('\n')
                    pre_exp = self.Reaction[i].variant_pre_expon[j]
                    if SDBool:
                        StiffCorrCounter += 1
                        if self.scaledown_factors[StiffCorrCounter] != 1:
                            txt.write('    {:25}{:.5e}'.format('pre_expon',
                                      pre_exp))
                            txt.write(('    # Pre-exponential has been ' +
                                      'rescaled by a factor of {0:.5e}\n').
                                      format(self.scaledown_factors
                                             [StiffCorrCounter]))
                        else:
                            txt.write('    {:25}{:.5e}\n'.format('pre_expon',
                                      pre_exp))
                    else:
                        txt.write('    {:25}{:.5e}\n'.format('pre_expon',
                                                             pre_exp))
                    if self.Reaction[i].reversible:
                        txt.write('    {:25}{:.5e}\n'.format('pe_ratio',
                                  self.Reaction[i].variant_pe_ratio[j]))
                    txt.write('    {:25}{:4.2f}\n'.format('activ_eng',
                              (self.Reaction[i].activ_eng[j])))
                    if self.Reaction[i].variant_prox_factor != []:
                        txt.write('    {:25}{:5.3f}\n'.format('prox_factor',
                                  (self.Reaction[i].
                                      variant_prox_factor[j])))
                    txt.write('  end_variant\n\n')

                if self.Reaction[i].reversible:
                    txt.write('end_reversible_step\n\n')
                else:
                    txt.write('end_step\n\n')

            txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_mechanism')

    def ReadSimIn(self):
        '''
        Read simulation_input.dat
        '''
        def StateInc(i):
            if _re.search('off', i):
                state = 'off'
                inc = ''
            elif _re.search('on time', i):
                state = 'time'
                inc = _np.float(i.split()[3])
            elif _re.search('on event', i):
                state = 'event'
                inc = _np.int(i.split()[3])
            return (state, inc)

        with open(_os.path.join(self.Path,
                                'simulation_input.dat'), 'r') as txt:
            RawTxt = txt.readlines()

        self.Conditions = _ClassIn.SimIn()
        self.Conditions.restart = True
        for i in RawTxt:
            if len(i.split()) > 0:
                if i[0] != '#':
                    i = i.split('#')[0]  # Don't parse comments
                    if i.split()[0] == 'temperature':
                        self.Conditions.T = _np.float(i.split()[1])
                    elif i.split()[0] == 'pressure':
                        self.Conditions.P = _np.float(i.split()[1])
                    elif i.split()[0] == 'random_seed':
                        self.Conditions.Seed = _np.int(i.split()[1])
                    elif i.split()[0] == 'no_restart':
                        self.Conditions.restart = False
                    elif i.split()[0] == 'gas_specs_names':
                        self.Conditions.gas_spec = i.split()[1:]
                        self.Conditions.n_gas = len(self.Conditions.gas_spec)
                    elif i.split()[0] == 'gas_energies':
                        for j in i.split()[1:]:
                            self.Conditions.gas_eng.append(_np.float(j))
                    elif i.split()[0] == 'gas_molec_weights':
                        for j in i.split()[1:]:
                            self.Conditions.gas_MW.append(_np.float(j))
                    elif i.split()[0] == 'gas_molar_fracs':
                        for j in i.split()[1:]:
                            self.Conditions.gas_molfrac.append(_np.float(j))
                    elif i.split()[0] == 'surf_specs_names':
                        self.Conditions.surf_spec = i.split()[1:]
                        self.Conditions.n_surf = len(self.Conditions.surf_spec)
                    elif i.split()[0] == 'surf_specs_dent':
                        for j in i.split()[1:]:
                            self.Conditions.surf_dent.append(_np.int(j))
                    elif i.split()[0] == 'event_report':
                        self.Conditions.event = i.split()[1]
                    elif i.split()[0] == 'snapshots':
                        self.Conditions.hist = StateInc(i)
                    elif i.split()[0] == 'process_statistics':
                        self.Conditions.procstat = StateInc(i)
                    elif i.split()[0] == 'species_numbers':
                        self.Conditions.specnum = StateInc(i)
                    elif i.split()[0] == 'max_time':
                        if i.split()[1] == 'infinity':
                            self.Conditions.SimTime_Max = 'inf'
                        else:
                            self.Conditions.SimTime_Max =\
                             _np.float(i.split()[1])
                    elif i.split()[0] == 'max_steps':
                        if i.split()[1] == 'infinity':
                            self.Conditions.MaxStep = 'inf'
                        else:
                            self.Conditions.MaxStep = int(i.split()[1])
                    elif i.split()[0] == 'wall_time':
                        self.Conditions.WallTime_Max = _np.int(i.split()[1])
                    elif i.split()[0] == 'finish' or\
                                         i.split()[0] == 'n_gas_species' or\
                                         i.split()[0] == 'n_surf_species':
                        pass

    def WriteSimIn(self):
        '''
        Write simulation_input.dat
        '''
        if not hasattr(self, 'Conditions'):
            self.ReadSimIn()

        with open(_os.path.join(self.Path,
                                'simulation_input.dat'), 'w') as txt:
            SeedTxt = ''
            if self.Conditions.Seed is None:
                _random.seed()
                self.Conditions.Seed = _random.randint(10000, 99999)
                SeedTxt = '# Random seed from Python wrapper'

            txt.write('#KMC simulation specification\n\n')
            txt.write('{:20}{:15}{}\n\n'.format('random_seed',
                      str(self.Conditions.Seed), SeedTxt))
            txt.write('{:20}{:5g}\n'.format('temperature', self.Conditions.T))
            txt.write('{:20}{:5g}\n\n'.format('pressure', self.Conditions.P))
            txt.write('{:20}{}\n'.format('n_gas_species',
                      str(self.Conditions.n_gas)))
            txt.write('{:20}'.format('gas_specs_names'))
            for i in range(0, self.Conditions.n_gas):
                txt.write('{:15}'.format(self.Conditions.gas_spec[i]))

            GasList = ['gas_energies', 'gas_molec_weights', 'gas_molar_fracs']
            GasList2 = ['gas_eng', 'gas_MW', 'gas_molfrac']
            for j in range(0, len(GasList)):
                txt.write('\n{:20}'.format(GasList[j]))
                for i in range(0, self.Conditions.n_gas):
                    txt.write('{:15}'.format(str(getattr(self.Conditions,
                              GasList2[j])[i])))

            txt.write('\n\n{:20}{}\n'.format('n_surf_species',
                      str(self.Conditions.n_surf)))
            txt.write('{:20}'.format('surf_specs_names'))
            for i in range(0, self.Conditions.n_surf):
                txt.write('{:15}'.format(self.Conditions.surf_spec[i]))
            txt.write('\n{:20}'.format('surf_specs_dent'))

            for i in range(0, self.Conditions.n_surf):
                txt.write('{:15}'.format(str(self.Conditions.surf_dent[i])))
            txt.write('\n\n')

            if self.Conditions.hist[0] == 'off':
                txt.write('{:20}{}\n'.format('snapshots', 'off'))
            elif self.Conditions.hist[0] == 'event':
                txt.write('{:20}{} {} {}\n'.format('snapshots', 'on',
                          self.Conditions.hist[0],
                          str(int(self.Conditions.hist[1]))))
            elif self.Conditions.hist[0] == 'time':
                txt.write('{:20}{} {} {}\n'.format('snapshots', 'on',
                          self.Conditions.hist[0],
                          str(_np.float(self.Conditions.hist[1]))))
            if self.Conditions.procstat[0] == 'off':
                txt.write('process_statistics  off\n')
            elif self.Conditions.procstat[0] == 'event':
                txt.write('{:20}{} {} {}\n'.format('process_statistics', 'on',
                          self.Conditions.procstat[0],
                          str(int(self.Conditions.procstat[1]))))
            elif self.Conditions.procstat[0] == 'time':
                txt.write('{:20}{} {} {}\n'.format('process_statistics', 'on',
                          self.Conditions.procstat[0],
                          str(_np.float(self.Conditions.procstat[1]))))

            if self.Conditions.specnum[0] == 'off':
                txt.write('species_numbers     off\n')
            elif self.Conditions.specnum[0] == 'event':
                txt.write('{:20}{} {} {}\n'.format('species_numbers', 'on',
                          self.Conditions.specnum[0],
                          str(int(self.Conditions.specnum[1]))))
            elif self.Conditions.specnum[0] == 'time':
                txt.write('{:20}{} {} {}\n'.format('species_numbers', 'on',
                          self.Conditions.specnum[0],
                          str(_np.float(self.Conditions.specnum[1]))))
            txt.write('{:20}{}\n\n'.format('event_report',
                      self.Conditions.event))

            if self.Conditions.MaxStep == '' or\
               _re.search('inf', str(self.Conditions.MaxStep)):
                txt.write('{:20}{}\n'.format('max_steps', 'infinity'))
            else:
                txt.write('{:20}{}\n'.format('max_steps',
                          str(self.Conditions.MaxStep)))

            if self.Conditions.SimTime_Max == '' or\
               _re.search('inf', str(self.Conditions.SimTime_Max)):
                txt.write('{:20}{}\n'.format('max_time', 'infinity\n'))
            else:
                txt.write('{:20}{}\n'.format('max_time',
                          str(self.Conditions.SimTime_Max)))
            if self.Conditions.WallTime_Max == '' or\
               _re.search('inf', str(self.Conditions.WallTime_Max)):
                txt.write('\n')
            else:
                txt.write('\n{:20}{}\n\n'.format('wall_time',
                          str(self.Conditions.WallTime_Max)))

            if not self.Conditions.restart:
                txt.write('no_restart\n')
            txt.write('finish\n')

    def ReadLatticeIn(self):
        '''
        Read lattice_input.dat
        '''
        self.Lattice = _ClassIn.LatticeIn()
        with open(_os.path.join(self.Path,
                                'lattice_input.dat'), 'r') as Txt:
            RawTxt = Txt.readlines()
        for i in RawTxt:
            self.Lattice.input.append(i.split('\n')[0])

    def WriteLattice(self):
        '''
        Write lattice_input.dat
        '''
        if not hasattr(self, 'Lattice'):
            self.ReadLatticeIn()

        with open(_os.path.join(self.Path,
                                'lattice_input.dat'), 'w') as txt:
            for i in self.Lattice.input:
                txt.write(i + '\n')

    def ReadStateInput(self):
        '''
        Read state_input.dat
        '''
        self.State = _ClassIn.StateIn()
        with open(_os.path.join(self.Path,
                                'state_input.dat'), 'r') as Txt:
            RawTxt = Txt.readlines()
        for i in RawTxt:
            self.State.Struct.append(i.split('\n')[0])
        self.State.Type = 'StateInput'

    def WriteStateIn(self):
        '''
        Write state_input.dat
        '''
        if not hasattr(self, 'State'):
            self.ReadStateInput()
        if self.State.Type == 'none':
            return

        if self.State.Type == 'StateInput':
            with open(_os.path.join(self.Path,
                                    'state_input.dat'), 'w') as txt:
                for i in self.State.Struct:
                    txt.write(i + '\n')

        elif self.State.Type == 'history':

            self.Lattice = self.State.Struct
            UniqSpec = _np.unique(self.Lattice[_np.not_equal(
                    self.Lattice[:, 2], 0), 1])
            nAds = len(UniqSpec)
            SpecIden = [0] * nAds
            AdsInfo = [[] for i in range(0, nAds)]
            DentInfo = [[] for i in range(0, nAds)]
            for i in range(0, nAds):
                for j in range(0, self.Lattice.shape[0]):
                    if UniqSpec[i] == self.Lattice[j, 1]:
                        AdsInfo[i].append(j + 1)
                        DentInfo[i].append(self.Lattice[j, 3])
                        SpecIden[i] = self.Lattice[j, 2]

            if nAds > 0:
                with open(_os.path.join(self.Path,
                                        'state_input.dat'), 'w') as txt:
                    txt.write('initial_state\n')
                    for i in range(0, nAds):
                        txt.write('  seed_on_sites  {:10}'.
                                  format(self.Conditions.surf_spec
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

    def ReadGeneral(self):
        '''
        Read general_output.txt
        '''
        with open(_os.path.join(self.Path, 'general_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()

        self.Performance = PerformanceIn()
        for i in range(0, len(RawTxt)):
            if _re.search('Number of elementary steps:', RawTxt[i]):
                nRxn = _np.int(RawTxt[i].split(':')[1])
            elif _re.search('Current KMC time:', RawTxt[i]):
                self.Performance.t_final = _np.float(RawTxt[i].split(':')[1])
            elif _re.search('Events occurred:', RawTxt[i]):
                self.Performance.events_occurred =\
                 _np.float(RawTxt[i].split(':')[1])
            elif _re.search('Elapsed CPU time:', RawTxt[i]):
                after_colon = RawTxt[i].split(':')[1]
                self.Performance.CPU_time =\
                    _np.float(after_colon.split(' ')[-2])
            elif _re.search('Reaction network:', RawTxt[i]):
                RxnStartLine = i + 2

        if RawTxt[RxnStartLine].split()[0] == '1.':
            NameInd = 1
        else:
            NameInd = 0

        RxnNameList = []
        nuList = []
        for i in range(RxnStartLine, RxnStartLine + nRxn):
            RxnName = RawTxt[i].split()[NameInd][:-1]
            RxnNameList.append(RxnName)
            RxnStr = RawTxt[i][_re.search('Reaction:', RawTxt[i]).end():]
            RxnStrList = RxnStr.split()
            nu = [0] * (self.Conditions.n_surf + self.Conditions.n_gas)
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
                                   range(0, len(self.Conditions.surf_spec))
                                   if SurfIden ==
                                   self.Conditions.surf_spec[k]][0]
                        nu[SurfInd] += Sign
                elif RxnStrList[j] != '->' and RxnStrList[j] != '+':
                    GasInd = [k for k in
                              range(0, len(self.Conditions.gas_spec))
                              if RxnStrList[j] ==
                              self.Conditions.gas_spec[k]][0]
                    nu[self.Conditions.n_surf + GasInd] += Sign
            nuList.append(nu)

        self.Performance.Nu = nuList
        self.Performance.UniqNu = _ReturnUnique(nuList).tolist()

    def ReadHistory(self):
        '''
        Read history_output.txt
        '''
        # Check if file exists
        if not _os.path.isfile(_os.path.join(self.Path, 'history_output.txt')):
            return

        # Count number of sites
        with open(_os.path.join(self.Path, 'lattice_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()
        nSites = len(RawTxt) - 2
        self.n_sites = nSites
        HistPath = _os.path.join(self.Path, 'history_output.txt')
        nLines = _rawbigcount(HistPath)
        self.n_snapshots = (nLines-6)/(nSites+2)
        self.History = []
        self.snap_times = []

        for snap_ind in range(self.n_snapshots):
            snap_data = _np.array([[0]*4]*nSites)
            _linecache.clearcache()
            snap_header = _linecache.getline(HistPath, 8 + snap_ind *
                                             (nSites+2)-1).split()
            self.snap_times.append(snap_header[3])
            for i in range(0, nSites):
                snap_data[i, :] = _linecache.getline(HistPath, 8 + snap_ind *
                                                     (nSites+2)+i).split()
            self.History.append(snap_data)

    def ReadProcstat(self):
        '''
        Read procstat_output.txt
        '''
        MaxLen = _np.int(2e4)
        with open(_os.path.join(self.Path, 'procstat_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen * 3:  # Procstat uses 3 lines per outputs
            Spacing = _np.int(_np.floor((len(RawTxt) - 1)/(MaxLen*3)))
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
            t.append(_np.float(RawTxt2[i*3].split()[3]))
            eventsTemp = RawTxt2[i*3+2].split()[1:]
            for j in range(0, len(eventsTemp)):
                eventsTemp[j] = _np.int(eventsTemp[j])
            events.append(eventsTemp)

        self.Procstat = _ClassIn.ProcstatIn()
        self.Procstat.Spacing = Spacing
        self.Procstat.t = _np.asarray(t)
        self.Procstat.events = _np.asarray(events)

    def ReadSpecnum(self):
        '''
        Read specnum_output.txt
        '''
        MaxLen = _np.int(2e4)
        with open(_os.path.join(self.Path, 'specnum_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen:
            Spacing = _np.int(_np.floor((len(RawTxt)-1)/MaxLen))
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
            nEvents.append(_np.int(LineSplit[1]))
            t.append(_np.float(LineSplit[2]))
            T.append(_np.float(LineSplit[3]))
            E.append(_np.float(LineSplit[4]))
            specTemp = LineSplit[5:]
            for j in range(0, len(specTemp)):
                specTemp[j] = _np.int(specTemp[j])
            spec.append(specTemp)
        self.Specnum = _ClassIn.SpecnumIn()
        self.Specnum.Spacing = Spacing
        self.Specnum.nEvents = _np.asarray(nEvents)
        self.Specnum.t = _np.asarray(t)
        self.Specnum.T = _np.asarray(T)
        self.Specnum.E = _np.asarray(E)
        self.Specnum.spec = _np.asarray(spec)

    def ReadCluster(self):
        '''
        Read clusterocc.bin
        '''
        dt = _np.dtype(_np.int32)
        virtual_arr = _np.memmap(_os.path.join(self.Path, 'clusterocc.bin'),
                                 dt, "r")
        nCluster = sum(len(s.variant_name) for s in self.Cluster)
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nCluster)
        if not hasattr(self, 'Binary'):
            self.Binary = _ClassIn.BinaryIn()
        self.Binary = _ClassIn.BinaryIn()
        self.Binary.cluster = _np.array(_np.reshape(virtual_arr,
                                                    [nNum/nCluster, nCluster])
                                        [::self.Specnum.Spacing])
        del virtual_arr

    def ReadProp(self, Mode):
        '''
        Mode = 0: Read Prop_output.bin
        Mode = 1: Read PropCounter_output.bin
        '''
        dt = _np.dtype(_np.float64)
        if Mode == 0:     # Instantaneous propensities
            FileName = 'Prop_output.bin'
        elif Mode == 1:   # Integral propensities
            FileName = 'PropCounter_output.bin'

        virtual_arr = _np.memmap(_os.path.join(self.Path,
                                               FileName), dt, "r")
        nRxn = len(self.Performance.Nu)
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]
        if not hasattr(self, 'Binary'):
            self.Binary = _ClassIn.BinaryIn()

        if Mode == 0:
            self.Binary.prop = _np.reshape(virtual_arr, [nNum/nRxn, nRxn])
            self.Binary.prop = _np.array(self.Binary.prop
                                         [::self.Procstat.Spacing])
        if Mode == 1:
            self.Binary.propCounter = _np.reshape(virtual_arr,
                                                  [nNum/nRxn, nRxn])
            self.Binary.propCounter = _np.array(self.Binary.propCounter
                                                [::self.Procstat.Spacing])
        del virtual_arr

    def ReadSA(self):
        '''
        Read SA_output.bin
        '''
        dt = _np.dtype(_np.float64)
        FileName = 'SA_output.bin'
        if _os.path.isfile(_os.path.join(self.Path, FileName)):
            virtual_arr = _np.memmap(_os.path.join(self.Path,
                                                   FileName), dt, "r")
            nRxn = len(self.Performance.Nu)
            nNum = virtual_arr.shape[0]
            nNum = nNum - (nNum % nRxn)
            virtual_arr = virtual_arr[:nNum]
            if not hasattr(self, 'Binary'):
                self.Binary = _ClassIn.BinaryIn()

            self.Binary.W_sen_anal = _np.reshape(virtual_arr,
                                                 [nNum/nRxn, nRxn])
            self.Binary.W_sen_anal = _np.array(self.Binary.W_sen_anal
                                               [::self.Specnum.Spacing])

            del virtual_arr
        else:
            print 'No sensitivity analysis output file'

    def CalcThermo(self, T):
        '''
        Calculate the forward activation energy, forward and reverse
        pre-exponential factors and the PE-ratio for each reaction described
        in Mechanism_input.dat using an input file with energies and
        vibrational frequencies for all species and transition states
        '''
        filepath = _os.path.join(self.Path, 'Zacros_Species_Energy.txt')
        [lines, dict] = _thermo.DFTFileRead(filepath)
        T_species = []
        for s in lines[3:]:
            T_species.append(_thermo.Reference(s.split('\t'), dict,
                                               filepath, T))
        '''
        Create list of transition state species
        '''
        TST = []

        for y in range(0, len(T_species)):
            if T_species[y].name.startswith('TST') or T_species[y].name == '*':
                TST.append([T_species[y].name, y])

        for x in range(0, len(self.Reaction)):
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
            for e in self.Reaction[x].initial:
                surf_species.append(e.split())
            surf_prod = []
            for e in self.Reaction[x].final:
                surf_prod.append(e.split())
            activ_eng = 0.0
            fwd_pre = 0.0

            if TST_index == -1:
                '''
                Case = No transition state energetics provided
                '''
                if hasattr(self.Reaction[x], 'gas_reacs_prods'):
                    MW_gas = next(e.MW for e in T_species
                                  if e.name == self.Reaction[x].
                                  gas_reacs_prods[0])
                    q_vib_gas = next(e.q_vib for e in T_species
                                     if e.name == self.Reaction[x].
                                     gas_reacs_prods[0])
                    q_rot_gas = next(e.q_rot for e in T_species
                                     if e.name == self.Reaction[x].
                                     gas_reacs_prods[0])
                    q_trans2D_gas = next(e.q_trans2D for e in T_species
                                         if e.name == self.Reaction[x].
                                         gas_reacs_prods[0])
                    for y in range(0, len(surf_prod)):
                        if surf_prod[y][1] != '*' and\
                                              int(surf_prod[y][2]) == 1:
                            q_vib_surf.append(next(e.q_vib for e in T_species
                                                   if e.name ==
                                                   surf_prod[y][1]))

                if hasattr(self.Reaction[x], 'gas_reacs_prods') and\
                   int(self.Reaction[x].gas_reacs_prods[1]) == -1:
                    '''
                    No transition state and a gas reactant
                    Non-activated adsorbtion
                    '''
                    fwd_pre = T_species[x].A_st /\
                        _np.sqrt(2*_np.pi * MW_gas * _c.kb1*T)\
                        * 1e5 * self.scaledown_factors[x]

                    rev_pre = q_vib_gas * q_rot_gas * q_trans2D_gas /\
                        _np.product(q_vib_surf) * _c.kb1 * T/_c.h1

                elif hasattr(self.Reaction[x], 'gas_reacs_prods') and\
                        int(self.Reaction[x].gas_reacs_prods[1]) == 1:
                    '''
                    No transition state and a gas product
                    Non-activated desorbtion
                    '''
                    rev_pre = T_species[x].A_st /\
                        _np.sqrt(2*_np.pi * MW_gas * _c.kb1*T)\
                        * 1e5 * self.scaledown_factors[x]

                    fwd_pre = q_vib_gas * q_rot_gas * q_trans2D_gas /\
                        _np.product(q_vib_surf) *\
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
                    T_species[TST_Slab].etotal
                q_vib_TST = next(e.q_vib for e in T_species
                                 if e.name ==
                                 T_species[TST_index].name)
                if hasattr(self.Reaction[x], 'gas_reacs_prods'):
                    q_vib_gas = next(e.q_vib for e in T_species
                                     if e.name == self.Reaction[x].
                                     gas_reacs_prods[0])
                    q_rot_gas = next(e.q_rot for e in T_species
                                     if e.name == self.Reaction[x].
                                     gas_reacs_prods[0])
                    q_trans2D_gas = next(e.q_trans2D for e in T_species
                                         if e.name == self.Reaction[x].
                                         gas_reacs_prods[0])
                    A_st = next(e.A_st for e in T_species
                                if e.name ==
                                self.Reaction[x].gas_reacs_prods[0])
                    MW_gas = next(e.MW for e in T_species
                                  if e.name ==
                                  self.Reaction[x].gas_reacs_prods[0])
                    for y in range(0, len(surf_prod)):
                        if surf_prod[y][1] != '*' and\
                              int(surf_prod[y][2]) == 1:
                                q_vib_surf.append(next(e.q_vib
                                                  for e in T_species
                                                  if e.name ==
                                                  surf_prod[y][1]))
                    Q_gas = q_vib_gas * q_rot_gas * q_trans2D_gas

                if hasattr(self.Reaction[x], 'gas_reacs_prods') and\
                   int(self.Reaction[x].gas_reacs_prods[1]) == -1:
                    '''
                    Transition state and a gas reactant
                    Activated adsorbtion
                    '''
                    activ_eng -=\
                        next(e.etotal for e in T_species
                             if e.name == self.Reaction[x].gas_reacs_prods[0])
                    fwd_pre = q_vib_TST/Q_gas * A_st /\
                        _np.sqrt(2*_np.pi*MW_gas*_c.kb1*T)*1e5
                    rev_pre = q_vib_TST/_np.product(q_vib_surf) *\
                        (_c.kb1*T/_c.h1)
                elif hasattr(self.Reaction[x], 'gas_reacs_prods') and\
                        int(self.Reaction[x].gas_reacs_prods[1]) == 1:
                    '''
                    Transition state and a gas product
                    Activated desorbtion
                    '''
                    q_vib_surf = []
                    for y in range(0, len(surf_species)):
                        if surf_species[y][1] != '*' and\
                              int(surf_species[y][2]) == 1:
                                q_vib_surf.append(next(e.q_vib
                                                  for e in T_species
                                                  if e.name ==
                                                  surf_species[y][1]))
                    rev_pre = q_vib_TST/Q_gas * A_st /\
                        _np.sqrt(2*_np.pi*MW_gas*_c.kb1*T)*1e5
                    fwd_pre = q_vib_TST/_np.product(q_vib_surf) *\
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
                                (next(e.etotal for e in T_species
                                 if e.name == surf_species[y][1]) -
                                 T_species[TST_Slab].etotal)
                            q_vib_surf.append(next(e.q_vib for e in T_species
                                              if e.name == surf_species[y][1]))
                    q_vib_reactants = _np.product(q_vib_surf)
                    fwd_pre = q_vib_TST/q_vib_reactants * (_c.kb1*T/_c.h1)

                    q_vib_prod = []
                    for y in range(0, len(surf_prod)):
                        if int(surf_prod[y][2]) == 1 and\
                          surf_prod[y][1] != '*':
                            q_vib_prod.append(next(e.q_vib for e in T_species
                                                   if e.name ==
                                                   surf_prod[y][1]))
                    q_vib_products = _np.product(q_vib_prod)
                    rev_pre = q_vib_TST/q_vib_products * (_c.kb1*T/_c.h1)
            self.Reaction[x].activ_eng[0] = max(activ_eng, 0.0)
            self.Reaction[x].variant_pre_expon[0] = fwd_pre
            self.Reaction[x].variant_pe_ratio[0] = fwd_pre/rev_pre


'''
 Class definitions for:

     Zacros input file      Class object
     -----------------      ------------
     energetics_input.dat   ClusterIn
     mechanism_input.dat    MechanismIn
     simulation_input.dat   SimIN
     lattice_input.dat      LatticeIn
     state_input.dat        StateIn
'''


class ClusterIn(IOdata):
    def __init__(self):
        self.variant_name = []
        self.variant_site_types = []
        self.variant_graph_multiplicity = []
        self.variant_eng = []
        self.nVariant = 0
        super(ClusterIn, self).__init__()


class MechanismIn(IOdata):
    def __init__(self):
        self.initial = []
        self.final = []
        self.variant_name = []
        self.variant_site_types = []
        self.variant_pre_expon = []
        self.variant_pe_ratio = []
        self.activ_eng = []
        self.variant_prox_factor = []
        self.nVariant = 0
        super(MechanismIn, self).__init__()


class SimIn(IOdata):
    def __init__(self):
        self.gas_eng = []
        self.gas_MW = []
        self.gas_molfrac = []
        self.surf_dent = []
        self.Seed = None
        self.WallTime_Max = ''
        super(SimIn, self).__init__()


class LatticeIn(IOdata):
    def __init__(self):
        self.input = []
        super(LatticeIn, self).__init__()


class StateIn(IOdata):
    def __init__(self):
        self.Struct = []
        super(StateIn, self).__init__()


class PerformanceIn(IOdata):
    def __init__(self):
        super(PerformanceIn, self).__init__()


class ProcstatIn(IOdata):
    def __init__(self):
        super(ProcstatIn, self).__init__()


class SpecnumIn(IOdata):
    def __init__(self):
        super(SpecnumIn, self).__init__()


class BinaryIn(IOdata):
    def __init__(self):
        super(BinaryIn, self).__init__()
