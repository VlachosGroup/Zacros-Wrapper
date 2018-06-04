# -*- coding: utf-8 -*-
'''
         -----------------------------------------------------
          Calculate thermodynamic data (S298, H298, and Cp(T)
          from ab initio DFT data (energies and frequencies)
              providng input thermodynamics files for
             KMC (Zacros) and MKM (Chemkin and Matlab)

                     Vlachos Research Group
                Chemical and Biomolecular Egineering
                      University of Delaware

                     Gerhard R Wittreich, P.E.
          -----------------------------------------------------

Created on Fri Mar 31 16:20:43 2017

@author: wittregr

Adopted from Matlab code written and modified by:

                     Vassili Vorotnikov
                            and
                         Geun Ho Gu

 This program contains the class objects used to read energy, vibration and
 molecular configuration data and determine the standard entropy and enthalpy
 and heat capacities at various temperatures.

'''
import numpy as _np
from numpy import pi
import os
import ase.io as _ase
import re
import scipy.interpolate as _sp
import datetime
from utils import constant as c


class Particle(object):

    Cp_Range = _np.linspace(100, 1500, 15)
    VibScalingFactor = 1         # Vibrational Scaling Factor

    def __init__(self, data, dict, Base_path, Tstp=298.15):

        '''
        Fill object with species name and associated thermodynamic data
        '''
        self.name = str(data[dict['name']])         # Species name
        if 'numh' in dict:
            if data[dict['numc']] != '':
                self.carbon = int(data[dict['numc']])       # No. of C atoms
            else:
                self.carbon = int(0)
            if data[dict['numh']] != '':
                self.hydrogen = int(data[dict['numh']])     # No. of H atoms
            else:
                self.hydrogen = int(0)
            if data[dict['numo']] != '':
                self.oxygen = int(data[dict['numo']])       # No. of O atoms
            else:
                self.oxygen = int(0)
            if data[dict['numn']] != '':
                self.nitrogen = int(data[dict['numn']])     # No. of N atoms
            else:
                self.nitrogen = int(0)
        if 'mw' in dict:
            self.MW = float(data[dict['mw']])/c.NA/1000
        else:
            self.MW = (self.carbon * c.MW_carbon +
                       self.hydrogen * c.MW_hydorgen +
                       self.oxygen * c.MW_oxygen +
                       self.nitrogen * c.MW_nitrogen)/c.NA/1000
        if not hasattr(self, 'Inertia'):
            self.totengpath = os.path.join(*re.split(r'\\|/',
                                           str(data[dict['totengpath']]).
                                           strip('.').strip('\\')))
        self.etotal = float(data[dict['etotal']])   # Total energy
        if data[dict['edisp']] == '':
            self.edisp = float(0)
        else:
            self.edisp = float(data[dict['edisp']])     # Dispersion energy
        self.numvibfreq = int(data[dict['numvibfreq']])  # No. vib frequencies
        if not hasattr(self, 'phase'):
            self.phase = None                       # Phase (G=gas, S=surface)
        self.vibfreq = []                           # Vibration frequencies
        for x in range(0, self.numvibfreq):
            self.vibfreq.append(Particle.VibScalingFactor *
                                float(data[dict['vibfreq'] + x]))
        self.vibfreq = _np.array(self.vibfreq)
        self.Base_path = Base_path
        self.Tstp = Tstp
        self.ThermoProperties()

    def ThermoProperties(self):
        '''
        Calculate all thermodynamic properties from input data
        '''

        self.Determine_Phase()

        '''
        Get rotational data from VASP CONTCAR for gas species using
        Atomic Simulation Environment (ASE) libraries for python
        '''
        if self.phase == 'G':
            if hasattr(self, 'Inertia'):
                self.I3 = self.Inertia
            else:
                filepath = os.path.join(self.Base_path,
                                        self.totengpath,
                                        'CONTCAR')
                VASP = _ase.read(filepath)
                self.I3 = VASP.get_moments_of_inertia()*c.A2_to_m2*c.amu_to_kg
                self.MW = sum(VASP.get_masses())/c.NA/1000.
            self.T_I = c.h1**2/(8*_np.pi**2*c.kb1)
        '''
        Calulcate common frequency data for vibrational components
        '''
        self.nu = self.vibfreq * 100 * c.c2
        self.theta = c.h1 * self.nu / c.kb1

        '''
        Call Entropy method to calculate standard state entropy
        '''
        self.Calc_Entropy()

        '''
        Call Heat Capacity method to calculate heat capacities at the
        temperature range specified by Cp_Range
        '''
        self.Calc_HeatCapacities()

        '''
        Call Enthalpy method to calculate standard state enthalpy
        '''
        self.Calc_Enthalpy()

    def Determine_Phase(self):
        '''
        Determine species phase if one is not provided
        '''

        if self.phase is not None:
            pass
        elif hasattr(self, 'surface'):
            self.phase = 'S'

        elif self.islinear is not None:
            self.phase = 'G'

        elif hasattr(self, 'sigma'):
            self.phase = 'G'

        else:
            '''
            This should proabbaly result in an error condition vs supplying
            'S' as a default value
            '''
            self.phase = 'S'

    def Calc_Entropy(self):
        '''
        Calculate the vibrational, rotational and translational components and
        total entropy of a species at standard conditions
        '''

        '''
        Calculate vibrational component of entropy for gas and surface species
        '''
        T = self.Tstp
        self.S_Tstp_vib = c.R1 * sum((self.theta/T) /
                                     (_np.exp(self.theta/T)-1) -
                                     _np.log(1 - _np.exp(-self.theta/T)))
        self.q_vib = _np.product(_np.divide(1, (1 - _np.exp(-self.theta/T))))

        '''
        Calculate rotational and translational components of entropy
        '''
        if self.phase == 'G':
            '''
            Gas phase calculation
            '''
            if self.islinear == 0:
                '''
                Non-linear species
                '''
                I = _np.product(self.I3)
                self.S_Tstp_rot = c.R1*(3./2. + 1./2. *
                                        _np.log(pi*T**3/self.T_I**3*I) -
                                        _np.log(self.sigma))
                self.q_rot = _np.sqrt(_np.pi*I)/self.sigma *\
                    (T/self.T_I)**(3./2.)
            else:
                '''
                Linear species
                '''
                I = _np.max(self.I3)
                self.S_Tstp_rot = c.R1*(1. + _np.log(T/self.T_I*I) -
                                        _np.log(self.sigma))
                self.q_rot = (T*I/self.sigma)/self.T_I
            p = 100000  # Presure of 1 atm or 100000 Pa
            self.S_Tstp_trans = c.R1*(5./2. + 3./2. *
                                      _np.log(2.*pi*self.MW/c.h1**2) +
                                      5./2.*_np.log(c.kb1*T) -
                                      _np.log(p))
            if hasattr(self, 'A_st'):
                self.q_trans2D = self.A_st * (2*pi*self.MW*c.kb1*T)/c.h1**2
        else:
            '''
            Surface phase calculation
            '''
            self.S_Tstp_rot = 0.
            self.S_Tstp_trans = 0.
            self.q_rot = 0.
        '''
        Sum of all contributions to entropy for total entropy
        '''
        self.S_Tstp = self.S_Tstp_vib + self.S_Tstp_rot + self.S_Tstp_trans

    def Calc_HeatCapacities(self):
        '''
        Calculate the vibrational, rotational and translational components and
        total heat capacity of a species for a range of temperatures
        '''
        '''
        Calculate vibrational contribution to heat capacity for temperature
        range specified in Cp_Range for linear and non-linear species
        '''
        zz = []
        for x in range(0, len(self.Cp_Range)):
            zz.append(_np.divide(self.theta, self.Cp_Range[x]))
        zz = _np.array(zz).T
        self.Cp_vib = sum(_np.divide(_np.multiply(_np.power(zz, 2),
                                                  _np.exp(-zz)),
                                     _np.power(1-_np.exp(-zz), 2)))
        if self.phase == 'G':
            '''
            Translational and rotational Gas phase calculation
            '''
            self.Cp_trans = _np.array([3./2.]*len(self.Cp_Range))
            if self.islinear == 0:
                '''
                Non-Linear species
                '''
                self.Cp_rot = _np.array([3./2.]*len(self.Cp_Range))
            else:
                '''
                Linear species
                '''
                self.Cp_rot = _np.array([1.]*len(self.Cp_Range))
        else:
            '''
            Surface species
            '''
            self.Cp_rot = _np.array([0.]*len(self.Cp_Range))
            self.Cp_trans = _np.array([0.]*len(self.Cp_Range))
        '''
        Sum all contribution to heat capacity for total heat capapcity
        '''
        self.Cp = c.R1*(self.Cp_trans + self.Cp_rot + self.Cp_vib + 1)

    def Calc_Enthalpy(self):
        T = self.Tstp
        '''
        Calculate zero-point energy
        '''
        self.zpe = sum(_np.multiply(c.h2, self.nu)/2)*c.NA/1000

        '''
        Calculate vibrational component of enthalpy
        '''
        self.E_Tstp_vib = c.kb2 *\
            sum(_np.divide(self.theta*_np.exp(-self.theta/T),
                           (1 - _np.exp(-self.theta/T)))) *\
            c.NA/1000
        '''
        Calculate translational and rotational component of enthalpy
        '''
        if self.phase == 'G':
            '''
            Gas phase calculation
            '''
            self.E_Tstp_trans = 3./2.*c.R1*T/1000
            if self.islinear == 0:
                '''
                Non-linear species
                '''
                self.E_Tstp_rot = 3./2.*c.R1*T/1000
            else:
                '''
                Linear species
                '''
                self.E_Tstp_rot = 1.*c.R1*T/1000
            self.pv_Tstp = 1.*c.R1*T/1000
        else:
            '''
            Surface phase calculation
            '''
            self.E_Tstp_trans = 0.
            self.E_Tstp_rot = 0.
            self.pv_Tstp = 0.
        '''
        Sum all contribution to enthalpy for total enthalpy
        '''
        self.dfth = self.etotal*c.R1/c.R2*c.NA/1000 + self.zpe +\
            self.E_Tstp_vib + self.E_Tstp_trans + self.E_Tstp_rot +\
            self.pv_Tstp


class Reference(Particle):
    '''
    SubClass object to add specific fields for reference species
    '''
    def __init__(self, data, dict, Base_path, Tstp=298.15):
        if data[dict['sigma']] != '':
            self.sigma = int(data[dict['sigma']])            # Sigma
        else:
            self.sigma = int(0)
        if data[dict['islinear']] != '':
            self.islinear = int(data[dict['islinear']])   # Is molecule linear?
        else:
            self.linear = int(-1)
        if 'hf298nist' in dict:
            self.hf298nist = float(data[dict['hf298nist']])  # NIST Std enthpy
        if 'inertia' in dict:
            if data[dict['inertia']] != '':
                self.Inertia = float(data[dict['inertia']])
            else:
                self.Inertia = float(0)
        self.phase = str.upper(data[dict['phase']])      # Phase
        if 'a_st' in dict and data[dict['a_st']] != '':
            self.A_st = float(data[dict['a_st']])
        super(Reference, self).__init__(data, dict, Base_path, Tstp=298.15)

    @staticmethod
    def BasisSet(RefSpecies):
        A = []
        b_nist = []
        b_dfth = []
        for x in range(0, len(RefSpecies)):
            A.append([RefSpecies[x].carbon, RefSpecies[x].hydrogen,
                      RefSpecies[x].oxygen, RefSpecies[x].nitrogen])
            b_nist.append([RefSpecies[x].hf298nist])
            b_dfth.append([RefSpecies[x].dfth])
        ref = _np.linalg.lstsq(A, _np.array(b_nist) - _np.array(b_dfth))[0]
        return(ref)


class Target(Particle):
    '''
    SubClass object to add specific fields for target surface species
    '''
    def __init__(self, data, dict, Base_path, Tstp=298.15):
        self.surface = str(data[dict['surface']])          # Surface
        self.functional = str(data[dict['functional']])    # Functional
        self.kpoints = str(data[dict['kpoints']])          # k-Points
        self.vibfreqpath = str(data[dict['vibfreqpath']])  # Unused
        self.phase = None                                  # Phase
        super(Target, self).__init__(data, dict, Base_path, Tstp=298.15)

    @staticmethod
    def ReferenceDFT(Species, Surface, Basis):
        for x in range(0, len(Species)):
            Molecule = _np.array([Species[x].carbon, Species[x].hydrogen,
                                  Species[x].oxygen, Species[x].nitrogen])
            if Species[x].phase == 'G':
                Species[x].hf_Tstp = (Species[x].dfth +
                                      _np.dot(Molecule, Basis))[0]
                if hasattr(Species[x], 'edisp'):
                    Species[x].convedisp = (Species[x].edisp *
                                            c.ev_atom_2_kcal_mol)
            else:
                Slab = next((y for y in Surface if y.name ==
                             Species[x].surface),
                            None)
                if Slab is None:
                    print 'Error'
                else:
                    Species[x].hf_Tstp = (Species[x].dfth +
                                          _np.dot(Molecule, Basis) -
                                          Slab.etotal *
                                          c.ev_atom_2_kcal_mol)[0]
                    if hasattr(Species[x], 'edisp'):
                        Species[x].convedisp = (Species[x].edisp *
                                                c.ev_atom_2_kcal_mol -
                                                Slab.edisp *
                                                c.ev_atom_2_kcal_mol)
        return(Species)

    @staticmethod
    def CreateThermdat(Species, Base_path, Output):
        '''
        Calculate the seven coefficients for the NASA polynomials for two
        temperature ranges and output the results in a Chemkin format thermdat
        file
        '''
        T_mid = 500
        Tstp = Species[0].Tstp

        def HS_NASA(T, a):
            '''
            7-coefficient NASA polynomials for enthalpy and entropy
            '''
            Enthalpy = a[0] + a[1]*T/2 + a[2]*T**2/3 + \
                a[3]*T**3/4 + a[4]*T**4/5
            Entropy = a[0]*_np.log(T) + a[1]*T + a[2]*T**2/2 + \
                a[3]*T**3/3 + a[4]*T**4/4
            return[Enthalpy, Entropy]

        for x in range(0, len(Species)):
            T_rng_low = _np.linspace(min(Species[x].Cp_Range), T_mid, 1600)
            T_rng_high = _np.linspace(T_mid, max(Species[x].Cp_Range), 4000)
            T_func = _sp.InterpolatedUnivariateSpline(Species[x].Cp_Range,
                                                      Species[x].Cp/c.R1, k=4)
            '''
            Fit coefficients A1-A5 to heat capacity data
            '''
            Species[x].a_low = _np.polyfit(T_rng_low,
                                           T_func(T_rng_low), 4)[::-1]
            Species[x].a_high = _np.polyfit(T_rng_high,
                                            T_func(T_rng_high), 4)[::-1]
            '''
            Correct A1 high temperature range coefficient to eliminate
            discontinuity between high and low temperature range polynomials
            '''
            Species[x].a_high[0] = Species[x].a_high[0] + \
                (_np.polyval(Species[x].a_low[::-1], T_mid) -
                 _np.polyval(Species[x].a_high[::-1], T_mid))

            '''
            Determine A6 coefficient for enthalpy calculations
            '''
            a6_high = (Species[x].hf_Tstp/c.R1/Tstp*1000 -
                       HS_NASA(Tstp, Species[x].a_high)[0])*Tstp
            a6_low = (Species[x].hf_Tstp/c.R1/Tstp*1000 -
                      HS_NASA(Tstp, Species[x].a_low)[0])*Tstp
            '''
            Correct A6 high temperature range coefficient to eliminate
            discontinuity between high and low temperature range polynomials
            '''
            a6_high_delta = (HS_NASA(T_mid, Species[x].a_low)[0] +
                             a6_low/T_mid) - \
                            (HS_NASA(T_mid,
                             Species[x].a_high)[0] + a6_high/T_mid)
            a6_high = a6_high + a6_high_delta * T_mid
            Species[x].a_high = _np.append(Species[x].a_high, a6_high)
            Species[x].a_low = _np.append(Species[x].a_low, a6_low)

            '''
            Determine A7 coefficient for entropy calculations
            '''
            a7_high = Species[x].S_Tstp/c.R1 - \
                HS_NASA(Tstp, Species[x].a_high)[1]
            a7_low = Species[x].S_Tstp/c.R1 - \
                HS_NASA(Tstp, Species[x].a_low)[1]

            '''
            Correct A7 high temperature range coefficient to eliminate
            discontinuity between high and low temperature range polynomials
            '''
            a7_high_delta = (HS_NASA(T_mid, Species[x].a_low)[1] +
                             a7_low) - (HS_NASA(T_mid,
                                        Species[x].a_high)[1] + a7_high)
            a7_high = a7_high + a7_high_delta
            Species[x].a_high = _np.append(Species[x].a_high, a7_high)
            Species[x].a_low = _np.append(Species[x].a_low, a7_low)

        '''
        Write the species name, seven NASA coefficients for both a high and
        a low temperature range and other data in the Chemkin thermdat
        file format
        '''
        if os.path.isdir(os.path.join(Base_path, Output)) is False:
            os.mkdir(os.path.join(Base_path, Output))
        filepath = os.path.join(Base_path, Output, 'thermdat')
        fid = open(filepath, 'w')
        fid.truncate()
        '''
        Write thermdat file header
        '''
        fid.write('THERMO ALL\n')
        for s in range(0, _np.size(Species)):
            '''
            Write header line for each species on line 1
            '''
            fid.write('%-16s' % (Species[s].name))
            fid.write('%-8s' % (datetime.date.today().strftime('%Y%m%d')))
            fid.write('%1s%-4i' % ('C', Species[s].carbon))
            fid.write('%1s%-4i' % ('O', Species[s].oxygen))
            fid.write('%1s%-4i' % ('H', Species[s].hydrogen))
            fid.write('%1s%-4i' % ('N', Species[s].nitrogen))
            if Species[s].name.find('(S)'):
                fid.write('S')
            else:
                fid.write('G')
            fid.write('%10.0f%10.0f%8.0f' % (min(Species[x].Cp_Range),
                                             max(Species[x].Cp_Range),
                                             T_mid))
            fid.write('%6s%1i\n' % ('', 1))
            '''
            Write first five NASA coefficients for
            low temperature range on line 2
            '''
            for x in range(0, 5):
                fid.write('%15E' % (Species[s].a_low[x]))
            fid.write('%4s%1i\n' % ('', 2))
            '''
            Write final two NASA coefficients for
            low temperature range on line 2
            '''
            for x in range(0, 2):
                fid.write('%15E' % (Species[s].a_low[x+4]))
            '''
            Write first three NASA coeficients for
            high temperature range on line 3
            '''
            for x in range(0, 3):
                fid.write('%15E' % (Species[s].a_high[x]))
            fid.write('%4s%1i\n' % ('', 3))
            '''
            Write final four NASA coefficients for
            high temperature range on line 4
            '''
            for x in range(0, 4):
                fid.write('%15E' % (Species[s].a_high[x+3]))
            fid.write('%19s%1i\n' % ('', 4))
        '''
        Write file footer and close the file
        '''
        fid.write('END\n')
        fid.close()

        return(Species)


class Surface:
    '''
    Class object to populate slab energies for surfaces
    '''
    def __init__(self, data, dict):
        self.name = str(data[dict['name']])          # Surface name
        self.etotal = float(data[dict['etotal']])    # Total energy-DFT
        self.edisp = float(data[dict['edisp']])      # Dispersion energy-DFT


def DFTFileRead(filepath):
    fid = open(filepath, 'r')
    file = fid.read()
    lines = file.splitlines()
    dict_array = lines[2].lower().split('\t')
    dict = {}
    for x in range(0, len(dict_array)):
        dict[dict_array[x]] = x
    return(lines, dict)
