# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 20:14:53 2016

@author: mpnun
"""

from const import const
import numpy as np

class Reaction:
    
    def __init__(self): 
        
        self.reactants = []         # list of reactant species
        self.products = []          # list of product species
        self.TS = []                # single transition state
      
    def K_rxn(self):
        Q_prod = 1.0
        for prod_spec in self.products:
            Q_prod = Q_prod * prod_spec.Q()
            
        Q_react = 1.0
        for react_spec in self.reactants:
            Q_react = Q_react * react_spec.Q()
            
        return Q_prod / Q_react
        
    def TS_BEP(omega,b):        # set TS information using BEP
        delE = 0
        Ea = omega * delE + b
        return Ea