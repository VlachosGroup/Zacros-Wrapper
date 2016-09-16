# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 14:51:25 2016

@author: mpnun
"""

class const:
    
    def __init__(self):
        self.h = 6.626e-34;                          # J*s
        self.kB = 1.3806488e-23;                     # J/K
        self.c = 3.0e10;                               # cm/s
        self.T_stp = 278.15;                         # K
        self.eVtoJ = 1.60217657e-19;                 # conversion factor 
        self.Na = 6.0221413e23;                      # Avogadros's number
        self.Jtokcal = 0.000239005736;               # conversion factor
        self.R = 1.986e-3;                           # kcal/K/mol
        self.Sden = 2.49081e-9;                      # mol / cm^2
        self.Ast = 1.0 / self.Sden / 10000 / self.Na;  # m^2
        self.Patm = 10^5;                            # Pa