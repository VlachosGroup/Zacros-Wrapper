# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:02:26 2016

@author: mpnun
"""

import os
from surf_spec import surf_spec
from Species import Species
import numpy as np

os.system('cls')

x = surf_spec()

x.vibs = np.array([100.0, 200.0, 300.0, 400.0])
print x.Q_vib()

#y = Species()
#print y.name