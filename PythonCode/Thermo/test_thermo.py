# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:02:26 2016

@author: mpnun
"""

import os
from Species import Species
import numpy as np

os.system('cls')

x = Species('CO',phase='surface')
x.vibs = np.array([100.0, 200.0, 300.0, 400.0])
x.calc_Q_vib()
print x.Q_vib