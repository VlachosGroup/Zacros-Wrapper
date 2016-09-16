# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 20:29:01 2016

@author: mpnun
"""

from const import const
import numpy as np

class ReactionNetwork:
    
    def __init__(self): 
        
        self.species = []
        self.reactions = []
        self.stoich_mat = np.array([])              # stoichiometric matrix
      
    def PlotEnrgDiagram(self,rxn_list):
        # rxn_list lists the reactions that happen and in what order
        print 'Plotting reaction energy diagram'