# -*- coding: utf-8 -*-
"""
Created on Sat May 20 15:00:49 2017

@author: wittregr
"""

import sys
import os

HomePath = os.path.expanduser('~')
sys.path.append(os.path.join(HomePath,'Documents','GitHub'))
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Zacros-Wrapper'))
sys.path.append(os.path.join(HomePath,'Documents','GitHub', 'Zacros-Wrapper','zacros_wrapper'))

import zacros_wrapper as zw


x = zw.kmc_traj()

x.Path = os.path.join(os.getcwd(), 'input')
x.ReadAllInput()



T = 650
filepath = os.path.join(x.Path, "Zacros_Species_Energy.xlsx")
x.mechin.CalcThermo(filepath,T)


x.Path = os.path.join(os.getcwd(), 'output')
x.WriteAllInput()
