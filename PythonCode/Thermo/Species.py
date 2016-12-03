# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 14:47:36 2016

@author: mpnun
"""

from const import const
import numpy as np

class Species(object):
    
    def __init__(self, name = '', phase = 'surface', E = 0, vibs = []): 
        
        self.phase = phase      # surface, non-linear gas, linear gas, TS
        self.name = name
        self.MW = 16.0               # species molecular weight (amu)
        self.elems = []             # elemental composition
        self.Q = 1.0
        self.Q_vib = 1.0
        self.Q_rot = 1.0
        self.Q_trans = 1.0
        self.E = E
        self.E_ZPE = 0
        
        vibs = np.array(vibs)
        self.vibs = vibs.astype(float)              # cm^-1             
    
    def calc_ZPE(self):
        self.E_ZPE = self.E + np.sum(const.h * const.c * self.vibs) / 2
    
    def calc_Q(self):
        self.calc_Q_vib()
        self.calc_Q_rot()
        self.calc_Q_trans()
        self.Q = self.Q_vib * self.Q_rot * self.Q_trans     
    
    def calc_Q_vib(self):
        x = const.h * const.c * self.vibs / (const.kB * const.T_stp)
        q_conts = np.exp(-x/2) / (1 - np.exp(-x))
        self.Q_vib = np.prod(q_conts)
        
    def calc_Q_rot(self):
        # will differ for linear and non-linear molecules
        self.Q_rot = 1.0
        
    def calc_Q_trans(self):
        # Have a 2D option in here
        self.Q_trans =  1.0