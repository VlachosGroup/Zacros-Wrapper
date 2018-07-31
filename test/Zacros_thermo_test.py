# -*- coding: utf-8 -*-
"""
Created on Sat May 20 15:00:49 2017

@author: wittregr
"""

import sys
import os
#sys.path.append('C:\\Users\\wittregr\\Documents\\Python Scripts')
sys.path.append('C:\\Users\\wittregr\\Documents\\GitHub\\Zacros-Wrapper\\zacros_wrapper')
import zacros_wrapper as zw
x = zw.kmc_traj()
#x.Path = os.path.join('C:\\Users\\wittregr\\Documents\\Python Scripts', 'Input')
x.Path = os.path.join(os.getcwd(), 'Input')
x.ReadAllInput()
print(x.mechin.rxn_list[0].variant_list[0].pe_ratio)
print(x.mechin.rxn_list[3].variant_list[0].activ_eng)
T = 650
filepath = os.path.join(x.Path, "Zacros_Species_Energy.xlsx")
x.mechin.CalcThermo(filepath,T)
print(x.mechin.rxn_list[0].variant_list[0].pe_ratio)
print(x.mechin.rxn_list[3].variant_list[0].activ_eng)
x.Path = os.path.join(os.getcwd(), 'Output')
x.WriteAllInput()
