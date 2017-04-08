# -*- coding: utf-8 -*-
"""
Created on Fri Apr 07 09:39:17 2017

@author: wittregr
"""

import DFT_to_Thermochemistry as _thermo
import os

Base_path = 'C:\Users\wittregr\Documents\Python Scripts'
os.chdir(Base_path)

'''
 Read reference species and calculate all thermodynamic quantities
'''
fid = open('Input\Reference_set_info.txt', 'r')
file = fid.read()
lines = file.splitlines()
dict_array = lines[2].lower().split('\t')
dict = {}
for x in range(0, len(dict_array)):
    dict[dict_array[x]] = x
T_ref = []
for s in lines[3:]:
    T_ref.append(_thermo.Reference(s.split('\t'), dict))

'''
 Read target species and calculate all thermodynamic quantities
'''
fid = open('Input\Tobe_Referenced.txt', 'r')
file = fid.read()
lines = file.splitlines()
dict_array = lines[2].lower().split('\t')
dict = {}
for x in range(0, len(dict_array)):
    dict[dict_array[x]] = x
T_target = []
for s in lines[3:]:
    T_target.append(_thermo.Target(s.split('\t'), dict))
