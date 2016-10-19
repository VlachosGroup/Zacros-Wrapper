# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 12:01:24 2016

@author: mpnun
"""

import sys
import os

sys.path.insert(0, '../KMCsim')
from Stats import Stats

os.system('cls')

# Test statistics

data1 = [3.1, 2.9, 2.99, 3.5, 2.2, 2.9]
data2 = [2.9, 2.7, 3.1, 3.3, 2.0, 2.8]

print Stats.diff_ci(data1, data2)