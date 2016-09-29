# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 20:14:53 2016

@author: mpnun
"""

from const import const
import numpy as np
from Species import Species

class Reaction:
    
    def __init__(self): 
        
        self.reactants = []         # list of reactant species
        self.products = []          # list of product species
        self.TS = Species(name = 'empty')                    # transition state
        self.delE = 0
        self.Ea_fwd = 0
        self.Ea_bwd = 0
        self.A_fwd = 1.0e12
        self.A_rat = 1.0
      
    def K_rxn(self):
        Q_prod = 1.0
        for prod_spec in self.products:
            Q_prod = Q_prod * prod_spec.Q()
            
        Q_react = 1.0
        for react_spec in self.reactants:
            Q_react = Q_react * react_spec.Q()
            
        return Q_prod / Q_react
    
    def calc_delE(self):        # set TS information using BEP
    
        E_react = 0
        Q_react = 1.0
        for spec in self.reactants:
            E_react += spec.E_ZPE
            Q_react = Q_react * spec.Q
            
        E_prod = 0
        Q_prod = 1.0
        for spec in self.products:
            E_prod += spec.E_ZPE
            Q_prod = Q_prod * spec.Q
        
        self.delE = E_prod - E_react
        self.A_rat = Q_prod / Q_react
        
        if self.TS.E_ZPE - E_react < 0 or self.TS.E_ZPE - E_react < self.delE:
            self.TS.name = 'empty'
        
        if self.TS.name == 'empty':
            self.A_fwd = const.kB* const.T_stp / const.h
            self.Ea_fwd = np.max([0, self.delE])
            self.Ea_bwd = np.max([0, -self.delE])
        else:
            self.A_fwd = const.kB* const.T_stp / const.h * self.TS.Q / Q_react
            self.Ea_fwd = np.max([self.TS.E_ZPE - E_react, 0, self.delE])
            self.Ea_bwd = np.max([self.TS.E_ZPE - E_prod, 0, -self.delE])
    
    def TS_BEP(self,omega,b):        # set TS information using BEP
        delE = 0
        Ea = omega * delE + b
        return Ea