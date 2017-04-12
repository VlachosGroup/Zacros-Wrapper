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


class Particle(object):

    h = 6.626070040E-34          # Planck constant [J-s]
    c = 299792458                # Speed of light [m/s]
    N_A = 6.022140857E23         # Avogadro's number [1/mol]
    kB = 1.38064852E-23          # Boltzmann constant [J/K]
    R = 1.9872036                # Gas constant [cal/mol-K]
    Tstp = 298.15                # Standard reference temperature [K]
    Cp_Range = _np.linspace(100, 1500, 15)
    VibScalingFactor = 1         # Vibrational Scaling Factor
    DFT_library_path = 'Input\Reference_set_info.txt'

    def __init__(self, data, dict, Base_path):

        '''
        Fill object with species name and associated thermodynamic data
        '''
        self.name = str(data[dict['name']])      # Species name [string]
        self.carbon = int(data[dict['numc']])    # No. of carbon atoms [int]
        self.hydrogen = int(data[dict['numh']])  # No. of hydrogem atoms [int]
        self.oxygen = int(data[dict['numo']])    # No. of oxygen atoms [int]
        self.nitrogen = int(data[dict['numn']])  # Number of nitrogen atoms
        self.totengpath = str(data[dict['totengpath']])  # VASP files location
        self.etotal = float(data[dict['etotal']])  # Total energy
        self.edisp = float(data[dict['edisp']])  # Dispersion energy
        self.numvibfreq = int(data[dict['numvibfreq']])  # No. vib frequencies
        self.phase = None                        # Phase (G=gas, S=surface)
        self.vibfreq = []                        # Vibration frequencies
        for x in range(0, self.numvibfreq):
            self.vibfreq.append(Particle.VibScalingFactor *
                                float(data[dict['vibfreq'] + x]))
        self.vibfreq = _np.array(self.vibfreq)
        self.Base_path = Base_path

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
            amu_to_kg = 1.66053904e-27  # kg/amu
            A2_to_m2 = 1.0e-20          # m^2/A^2
            self.totengpath = os.path.join(*re.split(r'\\|/', self.totengpath.strip('.').strip('\\')))
            filepath = os.path.join(self.Base_path.strip('.').strip('\\'),
                                    'Input', self.totengpath,
                                    'CONTCAR')
            VASP = _ase.read(filepath)
            self.I3 = VASP.get_moments_of_inertia()*A2_to_m2*amu_to_kg
            self.T_I = self.h**2/(8*pi**2*self.kB)
            self.MW = sum(VASP.get_masses())/self.N_A/1000.
        '''
        Calulcate common frequency data for vibrational components
        '''
        self.nu = self.vibfreq * 100 * self.c
        self.theta = self.h * self.nu / self.kB

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
        self.S_Tstp_vib = self.R * sum((self.theta/T)/(_np.exp(self.theta/T)-1) -
                                       _np.log(1 - _np.exp(-self.theta/T)))
        '''
        Calculate rotational and translational components of enthalpy
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
                self.S_Tstp_rot = self.R*(3./2. + 1./2.*_np.log(pi*T**3/self.T_I**3*I) -
                                          _np.log(self.sigma))
                self.Srotm = self.R*(1./2.*_np.log(_np.pi*T**3/self.T_I**3*I) -
                                     _np.log(self.sigma))
            else:
                '''
                Linear species
                '''
                I = _np.max(self.I3)
                self.S_Tstp_rot = self.R*(1. + _np.log(T/self.T_I*I) -
                                          _np.log(self.sigma))
                self.Srotm = self.R*(_np.log(T/self.T_I*I) -
                                     _np.log(self.sigma))
            p = 100000  # Presure of 1 atm or 100000 Pa
            self.S_Tstp_trans = self.R*(5./2. + 3./2.*_np.log(2.*pi*self.MW/self.h**2) +
                                        5./2.*_np.log(Particle.kB*T) - _np.log(p))
            self.Stransm = self.R*(3./2.*_np.log(2.*pi*self.MW/self.h**2) +
                                   5./2.*_np.log(Particle.kB*T) - _np.log(p))
        else:
            '''
            Surface phase calculation
            '''
            self.S_Tstp_rot = 0.
            self.S_Tstp_trans = 0.
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
        self.Cp = self.R*(self.Cp_trans + self.Cp_rot + self.Cp_vib + 1)

    def Calc_Enthalpy(self):
        h = 1.5836687E-34
        kB = 3.299829159077406E-24
        R2 = 5.189478952438986e+19  # eV/mol-K
        T = self.Tstp
        '''
        Calculate zero-point energy
        '''
        self.zpe = sum(_np.multiply(h, self.nu)/2)*self.N_A/1000
        '''
        Calculate vibrational component of enthalpy
        '''
        self.E_Tstp_vib = kB*sum(_np.divide(self.theta*_np.exp(-self.theta/T),
                                            (1 - _np.exp(-self.theta/T)))) *\
                                            self.N_A/1000
        '''
        Calculate translational and rotational component of enthalpy
        '''
        if self.phase == 'G':
            '''
            Gas phase calculation
            '''
            self.E_Tstp_trans = 3./2.*self.R*T/1000
            if self.islinear == 0:
                '''
                Non-linear species
                '''
                self.E_Tstp_rot = 3./2.*self.R*T/1000
            else:
                '''
                Linear species
                '''
                self.E_Tstp_rot = 1.*self.R*T/1000
            self.pv_Tstp = 1.*self.R*T/1000
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
        self.dfth = self.etotal*self.R/R2*self.N_A/1000 + self.zpe +\
            self.E_Tstp_vib + self.E_Tstp_trans + self.E_Tstp_rot +\
            self.pv_Tstp


class Reference(Particle):
    '''
    SubClass object to add specific fields for reference species
    '''
    def __init__(self, data, dict, Base_path):
        self.sigma = int(data[dict['sigma']])            # Sigma
        self.islinear = int(data[dict['islinear']])      # Is molecule linear?
        self.hf298nist = float(data[dict['hf298nist']])  # NIST Std enthalpy
        self.phase = str.upper(data[dict['phase']])  # Phase (G=gas, S=surface)
        super(Reference, self).__init__(data, dict, Base_path) # Call superclass


class Target(Particle):
    '''
    SubClass object to add specific fields for target surface species
    '''
    def __init__(self, data, dict, Base_path):
        self.surface = str(data[dict['surface']])          # Surface
        self.functional = str(data[dict['functional']])    # Functional
        self.kpoints = str(data[dict['kpoints']])          # k-Points
        self.vibfreqpath = str(data[dict['vibfreqpath']])  # Unused
        self.phase = None                            # Phase (G=gas, S=surface)
        super(Target, self).__init__(data, dict, Base_path)   # Call superclass
