# -*- coding: utf-8 -*-
"""
Created on Fri Dec 02 16:03:23 2016

@author: mpnun
"""

import os
import sys
import numpy as np

sys.path.append('C:/Users/mpnun/Dropbox/Github/ZacrosWrapper/PythonCode')
import KMCsim as zw
#from KMCsim.Stats import Stats

os.system('cls')

#x = np.array([ [0.9134, 0.6324, 0.0975], [0.2785, 0.5469, 0.9575], [0.9649, 0.1576, 0.9706], [0.9572, 0.4854, 0.8003], [0.1419, 0.4218, 0.9157] ])
#y = Stats.cov_mat_ci(np.transpose(x))
#print y



#data = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
#
#print zw.Stats.cov_ci( data, data )

R = np.random.random_sample([3,7])

x = zw.Stats.cov_mat_ci(R)

print x['cov_mat']
print x['ci_mat']