# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:44:51 2016

@author: mpnun
"""

# Inputs information for species and reactions in the water gas-shift network, then draws the reaction energy diagram

import os
import numpy as np
from ReactionNetwork import ReactionNetwork
from Species import Species
from Reaction import Reaction

os.system('cls')

WGS_net = ReactionNetwork()

''' Define Species '''

E_slab = -11448.22

WGS_net.species.append(Species(name = 'CO(g)', phase='gas', E = -588.9694, vibs = [2143]))      # 1
WGS_net.species.append(Species(name = 'CO2(g)', phase='gas', E = -1025.4623, vibs = [667, 1333, 2349]))      # 2
WGS_net.species.append(Species(name = 'O2(g)', phase='gas', E = -867.2017, vibs = [1580]))      # 3
WGS_net.species.append(Species(name = 'H2(g)', phase='gas', E = -31.6629, vibs = [4190]))      # 4
WGS_net.species.append(Species(name = 'H2O(g)', phase='gas', E = -467.4604, vibs = [1595, 3657, 37560]))      # 5
WGS_net.species.append(Species(name = 'CO*', phase='surface', E = -12039.09 - E_slab, vibs = [367, 374, 455, 2010]))      # 6
WGS_net.species.append(Species(name = 'O*', phase='surface', E = -11882.98 - E_slab , vibs = [366, 394, 467]))      # 7
WGS_net.species.append(Species(name = 'H2O*', phase='surface', E = -11916.04 - E_slab, vibs = [241, 454, 554, 1559, 3483, 3606]))      # 8
WGS_net.species.append(Species(name = 'OH*', phase='surface', E = -11899.15 - E_slab, vibs = [494, 855, 3576]))      # 9
WGS_net.species.append(Species(name = 'H*', phase='surface', E = -11464.54 - E_slab , vibs = [692, 701, 986]))      # 10
WGS_net.species.append(Species(name = 'HCO*', phase='surface', E = -12054.46 - E_slab, vibs = [235, 333, 377, 535, 718, 1129, 1284, 2907]))      # 11
WGS_net.species.append(Species(name = 'HCOO*', phase='surface', E = -12489.87 - E_slab, vibs = [287, 322, 336, 733, 904, 1264, 1288, 1535, 2995]))      # 12
WGS_net.species.append(Species(name = 'COOH*', phase='surface', E = -12490.26 - E_slab, vibs = [267, 308, 434, 572, 655, 1076, 1206, 1675, 3341]))      # 13

for spec in WGS_net.species:
    spec.calc_ZPE()
#    print spec.E_ZPE
    
''' Define reactions '''

WGS_net.AddRxn([1],[6])             # 1
WGS_net.AddRxn([3],[7,7])           # 2
WGS_net.AddRxn([10,10],[4])         # 3
WGS_net.AddRxn([2],[6,7])           # 4
WGS_net.AddRxn([5],[8])             # 5
WGS_net.AddRxn([8],[9,10])          # 6
WGS_net.AddRxn([9],[7,10])          # 7
WGS_net.AddRxn([9,9],[7,8])         # 8
WGS_net.AddRxn([6,9],[13])          # 9
WGS_net.AddRxn([13],[10,2])         # 10
WGS_net.AddRxn([13,7],[9,2])        # 11
WGS_net.AddRxn([13,9],[8,2])        # 12
WGS_net.AddRxn([10,6],[11])               # 13
WGS_net.AddRxn([11,7],[12])               # 14
WGS_net.AddRxn([12],[2,10])               # 15
WGS_net.AddRxn([12,7],[2,9])               # 16
WGS_net.AddRxn([12,9],[2,8])               # 17

''' Plot reaction energy diagram '''

rxn_order = [1,5,6,9,10,3]
WGS_net.PlotEnrgDiagram(rxn_order)