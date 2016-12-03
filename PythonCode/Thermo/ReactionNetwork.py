# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 20:29:01 2016

@author: mpnun
"""

from const import const
import numpy as np
from Reaction import Reaction
import matplotlib.pyplot as plt
import matplotlib as mat

class ReactionNetwork:
    
    def __init__(self): 
        
        self.species = []
        self.reactions = []
        self.stoich_mat = np.array([])              # stoichiometric matrix
      
    def PlotEnrgDiagram(self,rxn_list):
        
        # Initial state
        E = 0
        state_ind = [0,1]
        state_eng = [E,E]
        state = 2
        rxn_names = []
        rxn_x = []
        rxn_y = []        
        
        # Iterate through reactions in the pathway
        for i in rxn_list:
            
            TS_eng = E + self.reactions[i-1].Ea_fwd
            final_eng = E + self.reactions[i-1].delE
            
            # Put TS in if the reaction is activated
            if not self.reactions[i-1].TS.name is 'empty':
                
                state_ind.append(state)
                state_ind.append(state+1)
                state += 2
                state_eng.append(TS_eng)
                state_eng.append(TS_eng)
                
                rxn_names.append(self.reactions[i-1].text)
                rxn_x.append(state-1)
                rxn_y.append(state-1)
            
            state_ind.append(state)
            state_ind.append(state+1)
            state += 2
            state_eng.append(final_eng)
            state_eng.append(final_eng)
            
            E = final_eng
        
        # Plotting
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 16
        
        plt.figure()
        
        plt.plot(state_ind, state_eng, 'o-', markersize = 12)   
        
        for ind in range(len(rxn_names)):
            plt.text(rxn_x[ind], rxn_y[ind], rxn_names[ind], size = 24)
        
        plt.xticks(size=24)
        plt.yticks(size=24)
        plt.xlabel('reaction coordinate',size=30)
        plt.ylabel('energy (eV)',size=30)
#        plt.legend(rxn_labels,loc=1,prop={'size':20},frameon=False)
        plt.show()
        
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)        
        
    def AddRxn(self,rcnt_list,prod_list, name = ''):
        rxn = Reaction()
        rxn.text = name
        for rcnt in rcnt_list:
            rxn.reactants.append(self.species[rcnt-1])
        for prod in prod_list:
            rxn.products.append(self.species[prod-1])       
        rxn.calc_delE()
        self.reactions.append(rxn)
        
    def WriteRxnInfo(self):
        with open('rxn_parameters.txt', 'w') as txt:
            txt.write('Reaction name \t TS name \t Afwd A_fwd/A_rev \t Ea (eV) \t delta_E (eV) \n')
            for rxn in self.reactions:
                txt.write(rxn.text + '\t' + rxn.TS.name + '\t {0:.3E}'.format(rxn.A_fwd) + '\t {0:.3E}'.format(rxn.A_rat) + '\t {0:.3f}'.format(rxn.Ea_fwd) + '\t {0:.3f} \n'.format(rxn.delE))