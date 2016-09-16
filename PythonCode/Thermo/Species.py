# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 14:47:36 2016

@author: mpnun
"""

from const import const
import numpy as np

class Species(object):
    
    def __init__(self,name,phase): 
        
        self.phase = phase      # surface, non-linear gas, linear gas, TS
        self.name = name
        self.MW = 16.0               # species molecular weight (amu)
        self.elems = []             # elemental composition
        self.Q = 1.0
        self.Q_elect = 1.0
        self.Q_vib = 1.0
        self.Q_rot = 1.0
        self.Q_trans = 1.0
        self.E = 0
        self.E_ZPE = 0
        self.vibs = []              # cm^-1
        
    def Q(self):
        self.Q_vib() = self.Q_elect * self.Q_vib * self.Q_rot * self.Q_trans
    
    def Q_elect(self):
        self.Q_elect = np.exp(-self.E / (const().kB * const().T_stp))        
    
    def ZPE(self):
        self.E_ZPE = self.E + np.sum(const().h * const().c * self.vibs)
    
    def Q_vib(self):
        x = const().h * const().c * self.vibs / (const().kB * const().T_stp)
        q_conts = np.exp(-x/2) / (1 - np.exp(-x)); 
        self.Qvib = np.prod(q_conts)
        
    def Q_rot(self):
        # will differ for linear and 
        self.Q_rot = 1.0
        
    def Q_trans(self):
        # Have a 2D option in here
        self.Q_trans =  1.0