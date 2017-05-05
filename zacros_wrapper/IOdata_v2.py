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
import IOdata_v2 as _ClassIn

from Helper import ReadWithoutBlankLines as _ReadWithoutBlankLines


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

        if _os.path.isfile(_os.path.join(self.Path, 'Input',
                                         'state_input.dat')):
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

    def ReadEngIn(self):
        '''
        Read energetics_input.dat
        '''
        RawTxt = _ReadWithoutBlankLines(_os.path.join(self.Path, 'Input',
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

        with open(_os.path.join(self.Path, 'Output',
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
        RawTxt = _ReadWithoutBlankLines(_os.path.join(self.Path, 'Input',
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
    #                elif not InVariant and i not in StateLine:
    #                    print 'Unparsed line in mechanism input:'
    #                    print RawTxt[i]
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
                        self.Reaction[j].variant_activ_eng.\
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
        with open(_os.path.join(self.Path, 'Output',
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
                    txt.write('    {:25}\n'.format('site_types'))
                    for k in range(0, len(self.Reaction[i].
                                          variant_site_types[j])):
                        txt.write('{} '.format(
                                self.Reaction[i].variant_site_types[j][k]))
                    txt.write('\n')
                    pre_exp = self.Reaction[i].variant_pre_expon[j]
                    if SDBool:
                        StiffCorrCounter += 1
                        if self.scaledown_factors[StiffCorrCounter] != 1:
                            txt.write('    {:25}{:.5e}\n'.format('pre_expon',
                                      pre_exp))
                            txt.write('    # Pre-exponential has been rescaled\
                                      by a factor of {0:.5e}\n'.
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
                    txt.write('    {:25}{}\n'.format('activ_eng',
                              str(self.Reaction[i].variant_activ_eng[j])))
                    if hasattr(self.Reaction[i], 'variant_prox_factor'):
                        txt.write('    {:25}{}\n'.format('prox_factor',
                                  str(self.Reaction[i].
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

        with open(_os.path.join(self.Path, 'Input',
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

        with open(_os.path.join(self.Path, 'Output',
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
        with open(_os.path.join(self.Path, 'Input',
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

        with open(_os.path.join(self.Path, 'Output',
                                'lattice_input.dat'), 'w') as txt:
            for i in self.Lattice.input:
                txt.write(i + '\n')

    def ReadStateInput(self):
        '''
        Read state_input.dat
        '''
        self.State = _ClassIn.StateIn()
        with open(_os.path.join(self.Path, 'Input',
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
            with open(_os.path.join(self.Path, 'Output',
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
                with open(_os.path.join(self.Path, 'Output',
                                        'state_input.dat'), 'w') as txt:
                    txt.write('initial_state\n')
                    for i in range(0, nAds):
                        txt.write('  seed_on_sites  {:10}'.
                                  format(self.Species['surf_spec']
                                         [SpecIden[i]-1], 10))
                        for j in range(0,len(DentInfo[i])):
                            for k in range(0,len(DentInfo[i])):
                                if j + 1 == DentInfo[i][k]:
                                    txt.write(str(AdsInfo[i][k]) + '  ')
                        txt.write('\n')
                    txt.write('end_initial_state\n')
        else:
            print 'Unrecognized state_input type'
            print 'state_input not written'

    def CheckComplete(self):
        '''
        Check to see if a Zacros run has completed successfully
        '''
        Complete = False
        if _os.path.isfile(_os.path.join(self.Path, 'general_output.txt')):
            with open(_os.path.join(self.Path,
                                    'general_output.txt'), 'r') as txt:
                RawTxt = txt.readlines()
            for i in RawTxt:
                if _re.search('Normal termination', i):
                    Complete = True
        return Complete

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
        self.variant_activ_eng = []
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
        super(SimIn, self).__init__()


class LatticeIn(IOdata):
    def __init__(self):
        self.input = []
        super(LatticeIn, self).__init__()


class StateIn(IOdata):
    def __init__(self):
        self.Struct = []
        super(StateIn, self).__init__()
