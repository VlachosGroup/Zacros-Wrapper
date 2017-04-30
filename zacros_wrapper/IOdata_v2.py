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

from Helper import ReadWithoutBlankLines

Path = 'C:\Users\wittregr\Documents\Python Scripts'

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


class ClusterIn(object):
    def __init__(self):
        self.variant_name = []
        self.variant_site_types = []
        self.variant_graph_multiplicity = []
        self.variant_eng = []


class MechanismIn(object):
    def __init__(self):
        self.initial = []
        self.final = []
        self.variant_name = []
        self.variant_site_types = []
        self.variant_pre_expon = []
        self.variant_pe_ratio = []
        self.variant_activ_eng = []
        self.variant_prox_factor = []


class SimIn(object):
    def __init__(self):
        self.gas_eng = []
        self.gas_MW = []
        self.gas_molfrac = []
        self.surf_dent = []
        self.Seed = None


class LatticeIn(object):
    def __init__(self):
        self.input = []


class StateIn(object):
    def __init__(self):
        self.Struct = []


def ReadAllInput():
    '''
    Read all input files. state_input.dat will be read only if it is present.
    '''
    Conditions = ReadSimIn()
    LatticeIn = ReadLatticeIn()
    Cluster = ReadEngIn()
    [Reaction, SD] = ReadMechIn()

    if _os.path.isfile(_os.path.join(Path, 'Input', 'state_input.dat')):
        State = ReadStateInput()
    else:
        State = None
    return (Conditions, LatticeIn, Cluster, Reaction, SD, State)


def WriteAllInput():
    '''
    Write all Zacros input files
    '''
    WriteSimIn()
    WriteMechanism()
    WriteEnergetics()
    WriteStateIn()
    WriteLattice()


def ReadEngIn():
    '''
    Read energetics_input.dat
    '''
    RawTxt = ReadWithoutBlankLines(_os.path.join(Path, 'Input',
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
    Cluster = []
    for j in range(0, nCluster):
        Cluster.append(ClusterIn())
        Cluster[j].name = RawTxt[ClusterInd[j, 0]].split()[1]
        Count = 0
        for i in range(ClusterInd[j, 0] + 1, ClusterInd[j, 1]):
            if RawTxt[i].split()[0] == 'variant':
                Count += 1
            elif RawTxt[i].split()[0] == 'sites':
                Cluster[j].sites = int(RawTxt[i].split()[1])
            elif RawTxt[i].split()[0] == 'neighboring':
                Cluster[j].neighboring = RawTxt[i].split()[1:]
            elif RawTxt[i].split()[0] == 'lattice_state':
                Cluster[j].latstate = RawTxt[i + 1:i + 1 + Cluster[j].sites]
                for k in range(0, len(Cluster[j].latstate)):
                    Cluster[j].latstate[k] = Cluster[j].latstate[k].\
                                             split('\n')[0]

        nVariant = Count
        nClusterTotal += nVariant
        variantInd = _np.array([[0, 0]]*nVariant)
        Count = 0
        for i in range(ClusterInd[j, 0]+1, ClusterInd[j, 1]):
            if RawTxt[i].split()[0] == 'variant':
                variantInd[Count, 0] = i
            if RawTxt[i].split()[0] == 'end_variant':
                variantInd[Count, 1] = i
                Count += 1

        for k in range(0, nVariant):
            for i in range(variantInd[k, 0], variantInd[k, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    Cluster[j].variant_name.append(RawTxt[i].split()[1])
                elif RawTxt[i].split()[0] == 'site_types':
                    Cluster[j].variant_site_types.append(RawTxt[i].split()[1:])
                elif RawTxt[i].split()[0] == 'graph_multiplicity':
                    Cluster[j].variant_graph_multiplicity.append(RawTxt[i].
                                                                 split()[1])
                elif RawTxt[i].split()[0] == 'cluster_eng':
                    Cluster[j].variant_eng.append(float(RawTxt[i].split()[1]))
    return(Cluster)


def WriteEnergetics(Cluster):
    '''
    Write energetics_input.dat
    '''
    nCluster = len(Cluster)

    with open(_os.path.join(Path, 'Output',
                            'energetics_input.dat'), 'w') as txt:
        txt.write('energetics\n\n')
        for i in range(0, nCluster):
            txt.write('#'*80 + '\n\n')
            txt.write('cluster ' + Cluster[i].name + '\n\n')
            txt.write('  sites ' + str(Cluster[i].sites) + '\n')

            if hasattr(Cluster[i], 'neighboring'):
                txt.write('  neighboring')
                for j in range(0, len(Cluster[i].neighboring)):
                    txt.write(' ' + Cluster[i].neighboring[j])
                txt.write('\n')

            txt.write('  lattice_state\n')
            for j in range(0, int(Cluster[i].sites)):
                txt.write(Cluster[i].latstate[j] + '\n')

            nVariant = len(Cluster[i].variant_name)
            txt.write('\n')
            for j in range(0, nVariant):
                txt.write('  {} {}\n'.format('variant',
                          Cluster[i].variant_name[j]))
                txt.write('    {:25}'.format('site_types'))
                for k in range(0, len(Cluster[i].variant_site_types[j])):
                    txt.write(Cluster[i].variant_site_types[j][k] + ' ')
                txt.write('\n')
                if int(Cluster[i].variant_graph_multiplicity[j]) > 0:
                    txt.write('    {:25}{}\n'.format('graph_multiplicity',
                              str(Cluster[i].variant_graph_multiplicity[j])))
                txt.write('    {:25}{}\n'.format('cluster_eng',
                          str(Cluster[i].variant_eng[j])))
                txt.write('  end_variant\n\n')

            txt.write('end_cluster\n\n')
        txt.write('#'*80 + '\n\n')
        txt.write('\n\nend_energetics')


def ReadMechIn():
    '''
    Read mechanism_input.dat
    '''
    RawTxt = ReadWithoutBlankLines(_os.path.join(Path, 'Input',
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
        scaledown_factors = [_np.float(i) for i in RawTxt[StiffCorrLine+2].
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

    Reaction = []
    for j in range(0, nMech):

        Reaction.append(MechanismIn())
        Reaction[j].reversible = Rxn[j]

        # Count the variants

        Reaction[j].name = RawTxt[MechInd[j, 0]].split()[1]
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
                Reaction[j].gas_reacs_prods = RawTxt[i].split()[1:]
            elif RawTxt[i].split()[0] == 'sites':
                nSites = int(RawTxt[i].split()[1])
                Reaction[j].sites = nSites
            elif RawTxt[i].split()[0] == 'neighboring':
                neighbor = RawTxt[i].split()[1:]
                Reaction[j].neighboring = neighbor
            elif RawTxt[i].split()[0] == 'initial':
                LatState = RawTxt[i+1:i+1+nSites]
                for k in range(0, len(LatState)):
                        Reaction[j].initial.append(LatState[k].split('\n')[0])
                for k in range(0, nSites):
                    StateLine.append(i+1+k)
            elif RawTxt[i].split()[0] == 'final':
                LatState = RawTxt[i + 1:i + 1 + nSites]
                for k in range(0, len(LatState)):
                        Reaction[j].final.append(LatState[k].split('\n')[0])
                for k in range(0, nSites):
                    StateLine.append(i+1+k)
#                elif not InVariant and i not in StateLine:
#                    print 'Unparsed line in mechanism input:'
#                    print RawTxt[i]
        nVariant = Count
        variantInd = _np.array([[0, 0]]*nVariant)
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

        for k in range(0, nVariant):
            for i in range(variantInd[k, 0], variantInd[k, 1]):
                if RawTxt[i].split()[0] == 'variant':
                    Reaction[j].variant_name.append(RawTxt[i].split()[1])
                elif RawTxt[i].split()[0] == 'site_types':
                    Reaction[j].variant_site_types.\
                     append(RawTxt[i].split()[1:])
                elif RawTxt[i].split()[0] == 'pre_expon':
                    Reaction[j].variant_pre_expon.\
                     append(float(RawTxt[i].split()[1]))
                elif RawTxt[i].split()[0] == 'pe_ratio':
                    Reaction[j].variant_pe_ratio.\
                     append(float(RawTxt[i].split()[1]))
                elif RawTxt[i].split()[0] == 'activ_eng':
                    Reaction[j].variant_activ_eng.\
                     append(float(RawTxt[i].split()[1]))
                elif RawTxt[i].split()[0] == 'prox_factor':
                    Reaction[j].variant_prox_factor.\
                     append(float(RawTxt[i].split()[1]))
                elif RawTxt[i].split()[0] == '#':
                    pass
#                    else:
#                        print 'Unparsed line in mechanism variant:'
#                        print RawTxt[i]
    pass
    if StiffCorrLine == -1:
        scaledown_factors = _np.ones(sum(len(s.variant_name)
                                     for s in Reaction))
    return(Reaction, scaledown_factors)


def WriteMechanism(Reaction=None, scaledown_factors=None):
    '''
    Write mechanism_input.dat
    '''
    if Reaction is None or scaledown_factors is None:
        print 'None'
        [Reaction, scaledown_factors] = ReadMechIn()

    if scaledown_factors is None:
        SDBool = False
    else:
        SDBool = True
    nMech = len(Reaction)
    StiffCorrCounter = -1
    with open(_os.path.join(Path, 'Output',
                            'mechanism_input.dat'), 'w') as txt:
        txt.write('mechanism\n\n')
        if SDBool:
            txt.write('# Automated stiffness reconditioning employed\n')
            txt.write('# \n')
            txt.write('# SDF: ')
            for i in scaledown_factors:
                txt.write('{0:.5e} \t'.format(i))
            txt.write('\n\n')
        for i in range(0, nMech):

            txt.write('#'*80 + '\n\n')

            if Reaction[i].reversible:
                txt.write('reversible_step ' + Reaction[i].name + '\n')
            else:
                txt.write('step ' + Reaction[i].name + '\n')

            txt.write('  sites ' + str(Reaction[i].sites) + '\n')
            if hasattr(Reaction[i], 'neighboring'):
                txt.write('  neighboring')
                for j in range(0, len(Reaction[i].neighboring)):
                    txt.write(' ' + Reaction[i].neighboring[j])
                txt.write('\n')

            if hasattr(Reaction[i], 'gas_reacs_prods'):
                txt.write('  gas_reacs_prods ' +
                          str(Reaction[i].gas_reacs_prods[0]) +
                          ' ' + str(Reaction[i].gas_reacs_prods[1]) + '\n')

            txt.write('  initial\n')
            for j in range(0, int(Reaction[i].sites)):
                txt.write(Reaction[i].initial[j] + '\n')

            txt.write('  final\n')
            for j in range(0, int(Reaction[i].sites)):
                txt.write(Reaction[i].final[j] + '\n')

            nVariant = len(Reaction[i].variant_name)
            txt.write('\n')
            for j in range(0, nVariant):
                txt.write('  variant ' + Reaction[i].variant_name[j] + '\n')
                txt.write('    {:25}\n'.format('site_types'))
                for k in range(0, len(Reaction[i].variant_site_types[j])):
                    txt.write(Reaction[i].variant_site_types[j][k] + ' ')
                txt.write('\n')
                pre_exp = Reaction[i].variant_pre_expon[j]
                if SDBool:
                    StiffCorrCounter += 1
                    if scaledown_factors[StiffCorrCounter] != 1:
                        txt.write('    {:25}{:.5e}\n'.format('pre_expon',
                                  pre_exp))
                        txt.write('    # Pre-exponential has been rescaled by\
                                  a factor of ' + '{0:.5e}'.
                                  format(scaledown_factors[StiffCorrCounter]) +
                                  ' \n')
                    else:
                        txt.write('    {:25}{:.5e}\n'.format('pre_expon',
                                  pre_exp))
                else:
                    txt.write('    {:25}{:.5e}\n'.format('pre_expon', pre_exp))
                if Reaction[i].reversible:
                    txt.write('    {:25}{:.5e}\n'.format('pe_ratio',
                              Reaction[i].variant_pe_ratio[j]))
                txt.write('    {:25}{}\n'.format('activ_eng',
                          str(Reaction[i].variant_activ_eng[j])))
                if hasattr(Reaction[i], 'variant_prox_factor'):
                    txt.write('    {:25}{}\n'.format('prox_factor',
                              str(Reaction[i].variant_prox_factor[j])))
                txt.write('  end_variant\n\n')

            if Reaction[i].reversible:
                txt.write('end_reversible_step\n\n')
            else:
                txt.write('end_step\n\n')

        txt.write('#'*80 + '\n\n')
        txt.write('\n\nend_mechanism')


def ReadSimIn():
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

    with open(_os.path.join(Path, 'Input',
                            'simulation_input.dat'), 'r') as txt:
        RawTxt = txt.readlines()

    Conditions = SimIn()
    Conditions.restart = True
    for i in RawTxt:
        if len(i.split()) > 0:
            if i[0] != '#':
                i = i.split('#')[0]  # Don't parse comments
                if i.split()[0] == 'temperature':
                    Conditions.T = _np.float(i.split()[1])
                elif i.split()[0] == 'pressure':
                    Conditions.P = _np.float(i.split()[1])
                elif i.split()[0] == 'random_seed':
                    Conditions.Seed = _np.int(i.split()[1])
                elif i.split()[0] == 'no_restart':
                    Conditions.restart = False
                elif i.split()[0] == 'gas_specs_names':
                    Conditions.gas_spec = i.split()[1:]
                    Conditions.n_gas = len(Conditions.gas_spec)
                elif i.split()[0] == 'gas_energies':
                    for j in i.split()[1:]:
                        Conditions.gas_eng.append(_np.float(j))
                elif i.split()[0] == 'gas_molec_weights':
                    for j in i.split()[1:]:
                        Conditions.gas_MW.append(_np.float(j))
                elif i.split()[0] == 'gas_molar_fracs':
                    for j in i.split()[1:]:
                        Conditions.gas_molfrac.append(_np.float(j))
                elif i.split()[0] == 'surf_specs_names':
                    Conditions.surf_spec = i.split()[1:]
                    Conditions.n_surf = len(Conditions.surf_spec)
                elif i.split()[0] == 'surf_specs_dent':
                    for j in i.split()[1:]:
                        Conditions.surf_dent.append(_np.int(j))
                elif i.split()[0] == 'event_report':
                    Conditions.event = i.split()[1]
                elif i.split()[0] == 'snapshots':
                    Conditions.hist = StateInc(i)
                elif i.split()[0] == 'process_statistics':
                    Conditions.procstat = StateInc(i)
                elif i.split()[0] == 'species_numbers':
                    Conditions.specnum = StateInc(i)
                elif i.split()[0] == 'max_time':
                    if i.split()[1] == 'infinity':
                        Conditions.SimTime_Max = 'inf'
                    else:
                        Conditions.SimTime_Max = _np.float(i.split()[1])
                elif i.split()[0] == 'max_steps':
                    if i.split()[1] == 'infinity':
                        Conditions.MaxStep = 'inf'
                    else:
                        Conditions.MaxStep = int(i.split()[1])
                elif i.split()[0] == 'wall_time':
                    Conditions.WallTime_Max = _np.int(i.split()[1])
                elif i.split()[0] == 'finish' or\
                                     i.split()[0] == 'n_gas_species' or\
                                     i.split()[0] == 'n_surf_species':
                    pass
    return Conditions


def WriteSimIn(Conditions=None):
    '''
    Write simulation_input.dat
    '''
    if Conditions is None:
        Conditions = ReadSimIn()

    with open(_os.path.join(Path, 'Output',
                            'simulation_input.dat'), 'w') as txt:
        SeedTxt = ''
        if Conditions.Seed is None:
            _random.seed()
            Conditions.Seed = _random.randint(10000, 99999)
            SeedTxt = '# Random seed from Python wrapper'

        txt.write('#KMC simulation specification\n\n')
        txt.write('{:20}{:15}{}\n\n'.format('random_seed',
                  str(Conditions.Seed), SeedTxt))
        txt.write('{:20}{:5g}\n'.format('temperature', Conditions.T))
        txt.write('{:20}{:5g}\n\n'.format('pressure', Conditions.P))
        txt.write('{:20}{}\n'.format('n_gas_species', str(Conditions.n_gas)))
        txt.write('{:20}'.format('gas_specs_names'))
        for i in range(0, Conditions.n_gas):
            txt.write('{:15}'.format(Conditions.gas_spec[i]))

        GasList = ['gas_energies', 'gas_molec_weights', 'gas_molar_fracs']
        GasList2 = ['gas_eng', 'gas_MW', 'gas_molfrac']
        for j in range(0, len(GasList)):
            txt.write('\n{:20}'.format(GasList[j]))
            for i in range(0, Conditions.n_gas):
                txt.write('{:15}'.format(str(getattr(Conditions,
                          GasList2[j])[i])))

        txt.write('\n\n{:20}{}\n'.format('n_surf_species',
                  str(Conditions.n_surf)))
        txt.write('{:20}'.format('surf_specs_names'))
        for i in range(0, Conditions.n_surf):
            txt.write('{:15}'.format(Conditions.surf_spec[i]))
        txt.write('\n{:20}'.format('surf_specs_dent'))

        for i in range(0, Conditions.n_surf):
            txt.write('{:15}'.format(str(Conditions.surf_dent[i])))
        txt.write('\n\n')

        if Conditions.hist[0] == 'off':
            txt.write('{:20}{}\n'.format('snapshots', 'off'))
        elif Conditions.hist[0] == 'event':
            txt.write('{:20}{} {} {}\n'.format('snapshots', 'on',
                      Conditions.hist[0], str(int(Conditions.hist[1]))))
        elif Conditions.hist[0] == 'time':
            txt.write('{:20}{} {} {}\n'.format('snapshots', 'on',
                      Conditions.hist[0], str(_np.float(Conditions.hist[1]))))
        if Conditions.procstat[0] == 'off':
            txt.write('process_statistics  off\n')
        elif Conditions.procstat[0] == 'event':
            txt.write('{:20}{} {} {}\n'.format('process_statistics', 'on',
                      Conditions.procstat[0],
                      str(int(Conditions.procstat[1]))))
        elif Conditions.procstat[0] == 'time':
            txt.write('{:20}{} {} {}\n'.format('process_statistics', 'on',
                      Conditions.procstat[0],
                      str(_np.float(Conditions.procstat[1]))))

        if Conditions.specnum[0] == 'off':
            txt.write('species_numbers     off\n')
        elif Conditions.specnum[0] == 'event':
            txt.write('{:20}{} {} {}\n'.format('species_numbers', 'on',
                      Conditions.specnum[0], str(int(Conditions.specnum[1]))))
        elif Conditions.specnum[0] == 'time':
            txt.write('{:20}{} {} {}\n'.format('species_numbers', 'on',
                      Conditions.specnum[0],
                      str(_np.float(Conditions.specnum[1]))))
        txt.write('{:20}{}\n\n'.format('event_report', Conditions.event))

        if Conditions.MaxStep == '' or\
           _re.search('inf', str(Conditions.MaxStep)):
            txt.write('{:20}{}\n'.format('max_steps', 'infinity'))
        else:
            txt.write('{:20}{}\n'.format('max_steps', str(Conditions.MaxStep)))

        if Conditions.SimTime_Max == '' or\
           _re.search('inf', str(Conditions.SimTime_Max)):
            txt.write('{:20}{}\n'.format('max_time', 'infinity\n'))
        else:
            txt.write('{:20}{}\n'.format('max_time',
                      str(Conditions.SimTime_Max)))
        if Conditions.WallTime_Max == '' or\
           _re.search('inf', str(Conditions.WallTime_Max)):
            txt.write('\n')
        else:
            txt.write('\n{:20}{}\n\n'.format('wall_time',
                      str(Conditions.WallTime_Max)))

        if not Conditions.restart:
            txt.write('no_restart\n')
        txt.write('finish\n')


def ReadLatticeIn():
    '''
    Read lattice_input.dat
    '''
    Lattice = LatticeIn()
    with open(_os.path.join(Path, 'Input', 'lattice_input.dat'), 'r') as Txt:
        RawTxt = Txt.readlines()
    for i in RawTxt:
        Lattice.Input.append(i.split('\n')[0])
    return Lattice


def WriteLattice(Lattice=None):
    '''
    Write lattice_input.dat
    '''
    if Lattice is None:
        Lattice = ReadLatticeIn()

    with open(_os.path.join(Path, 'Output', 'lattice_input.dat'), 'w') as txt:
        for i in Lattice.Input:
            txt.write(i + '\n')


def ReadStateInput():
    '''
    Read state_input.dat
    '''
    State = StateIn()
    with open(_os.path.join(Path, 'Input', 'state_input.dat'), 'r') as Txt:
        RawTxt = Txt.readlines()
    for i in RawTxt:
        State.Struct.append(i.split('\n')[0])
    State.Type = 'StateInput'
    return State


def WriteStateIn(State=None):
    '''
    Write state_input.dat
    '''
    if State is None:
        State = ReadStateInput()
    if State.Type == 'none':
        return

    if State.Type == 'StateInput':
        with open(_os.path.join(Path, 'Output',
                                'state_input.dat'), 'w') as txt:
            for i in State.Struct:
                txt.write(i + '\n')

    elif State.Type == 'history':

        Lattice = State.Struct
        UniqSpec = _np.unique(Lattice[_np.not_equal(Lattice[:, 2], 0), 1])
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
            with open(_os.path.join(Path, 'Output',
                                    'state_input.dat'), 'w') as txt:
                txt.write('initial_state\n')
                for i in range(0, nAds):
                    txt.write('  seed_on_sites  {:10}'.format(Species['surf_spec'][SpecIden[i]-1],10))
                    for j in range(0,len(DentInfo[i])):
                        for k in range(0,len(DentInfo[i])):
                            if j + 1 == DentInfo[i][k]:
                                txt.write(str(AdsInfo[i][k]) + '  ')
                    txt.write('\n')
                txt.write('end_initial_state\n')
    else:
        print 'Unrecognized state_input type'
        print 'state_input not written'
