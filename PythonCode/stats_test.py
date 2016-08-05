# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 22:53:42 2016

@author: mpnun
"""

from Stats import Stats

x = [1, 2, 3, 4, 5, 6]
y = [1.1, 1.9, 2.9, 3.5, 5.1, 8]

z = Stats().cov_calc(x,y)
print z
print Stats().cov_ci(x,y)