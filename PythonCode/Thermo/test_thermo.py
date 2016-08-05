# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:02:26 2016

@author: mpnun
"""

import os
from surf_spec import surf_spec
from Species import Species

os.system('cls')

x = surf_spec()

x.vibs = [100, 200, 300, 400]
print x.Q_vib()

#y = Species()
#print y.name