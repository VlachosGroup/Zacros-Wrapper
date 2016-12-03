# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:44:51 2016

@author: mpnun
"""

# Inputs information for species and reactions in the water gas-shift network, then draws the reaction energy diagram

import os
import sys
import numpy as np

sys.path.append('..')
from Thermo.ReactionNetwork import ReactionNetwork
from Thermo.Species import Species
from Thermo.Reaction import Reaction

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
    spec.calc_Q()
    
''' Define reactions '''


# Reactants and products

WGS_net.AddRxn([1],[6], name = 'CO(g) + * <-> CO*')             # 1
WGS_net.AddRxn([3],[7,7], name = 'O2(g) + 2 *<->2 O*')           # 2
WGS_net.AddRxn([10,10],[4], name = 'H2(g) + 2 *<->2 H*')         # 3
WGS_net.AddRxn([2],[6,7], name = 'CO2(g) + 2 *<->CO* + O*')           # 4
WGS_net.AddRxn([5],[8], name = 'H2O(g) + *<->H2O*')             # 5
WGS_net.AddRxn([8],[9,10], name = 'H2O* + *<->OH* + H*')          # 6
WGS_net.AddRxn([9],[7,10], name = 'OH* + *<->O* + H*')          # 7
WGS_net.AddRxn([9,9],[7,8], name = '2 OH*<->O* + H2O*')         # 8
WGS_net.AddRxn([6,9],[13], name = 'CO* + OH*<->COOH* +*')          # 9
WGS_net.AddRxn([13],[10,2], name = 'COOH*<->H* + CO2(g)')         # 10
WGS_net.AddRxn([13,7],[9,2], name = 'COOH* + O*<->OH* + CO2(g)')        # 11
WGS_net.AddRxn([13,9],[8,2], name = 'COOH* + OH*<->H2O* + CO2(g)')        # 12
WGS_net.AddRxn([10,6],[11], name = 'H* + CO*<->HCO**')               # 13
WGS_net.AddRxn([11,7],[12], name = 'HCO** + O*<->HCOO***')               # 14
WGS_net.AddRxn([12],[2,10], name = 'HCOO***<->CO2(g) + H* + 2 *')               # 15
WGS_net.AddRxn([12,7],[2,9], name = 'HCOO*** + O*<->CO2(g) + OH* + 3 *')               # 16
WGS_net.AddRxn([12,9],[2,8], name = 'HCOO*** + OH*<->CO2(g) + H2O* + 3 *')               # 17

# transitions states

WGS_net.reactions[2].TS = Species(name = 'TS3', phase='surface', E = -11479.95 - E_slab, vibs = [305, 310, 1601, 2255])
WGS_net.reactions[3].TS = Species(name = 'TS4', phase='surface', E = -12472.44 - E_slab, vibs = [290, 350, 370, 422, 545, 1896])
WGS_net.reactions[5].TS = Species(name = 'TS6', phase='surface', E = -11915.26 - E_slab, vibs = [254, 282, 493, 524, 870, 1866, 3552])
WGS_net.reactions[6].TS = Species(name = 'TS7', phase='surface', E = -11898.21 - E_slab, vibs = [260, 410, 477, 1796])
WGS_net.reactions[7].TS = Species(name = 'TS8', phase='surface', E = -12350.53 - E_slab, vibs = [268, 432, 565, 578, 724, 839, 1353, 1454, 3551])
WGS_net.reactions[8].TS = Species(name = 'TS9', phase='surface', E = -12489.54 - E_slab, vibs = [310, 386, 409, 512, 553, 839, 1926, 3545])
WGS_net.reactions[9].TS = Species(name = 'TS10', phase='surface', E = -12489.40 - E_slab, vibs = [275, 319, 555, 576, 1157, 1804, 1958])
WGS_net.reactions[10].TS = Species(name = 'TS11', phase='surface', E = -12924.54 - E_slab, vibs = [235, 264, 324, 448, 568, 630, 1048, 1120, 1174, 1764, 2234])
WGS_net.reactions[11].TS = Species(name = 'TS12', phase='surface', E = -12941.08 - E_slab, vibs = [267, 308, 434, 494, 572, 655, 855, 1076, 1206, 1675, 3341, 3576])
WGS_net.reactions[12].TS = Species(name = 'TS13', phase='surface', E = -12054.10 - E_slab, vibs = [213, 288, 315, 624, 977, 1597, 2053])
WGS_net.reactions[13].TS = Species(name = 'TS14', phase='surface', E = -12488.19 - E_slab, vibs = [285, 305, 412, 475, 555, 818, 1057, 1488, 2975])
WGS_net.reactions[14].TS = Species(name = 'TS15', phase='surface', E = -12488.83 - E_slab, vibs = [207, 370, 608, 895, 1030, 1240, 1532, 1644])
WGS_net.reactions[15].TS = Species(name = 'TS16', phase='surface', E = -12922.88 - E_slab, vibs = [235, 264, 324, 448, 568, 630, 1048, 1120, 1174, 1764, 2234])
WGS_net.reactions[16].TS = Species(name = 'TS17', phase='surface', E = -12939.78 - E_slab, vibs = [201, 225, 375, 473, 652, 674, 963, 977, 1103, 1292, 1621, 1959, 3052])

# Compute reaction info
for rxn in WGS_net.reactions:
    rxn.TS.calc_Q()
    rxn.calc_delE()  
            
''' Plot reaction energy diagram '''

rxn_order = [1,5,6,9,10,3]
WGS_net.PlotEnrgDiagram(rxn_order)
WGS_net.WriteRxnInfo()