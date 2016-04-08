# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 10:52:29 2016

@author: robieta
"""

import os
import numpy as np

class MachineSpecifics:
    def __init__(self):
        pass

    def MaxOutputEntries(self):
        # This function determines the maximum length of the arrays of output
        # files. This allows robust treatment of large files. Entries will be 
        # skipped to report the maximum KMC time subject to this constraint.
        MaxLen = 5e3
        MaxLen = np.int(MaxLen)
        return MaxLen
        
    def MinACFCoarseGrain(self):
        # This function sets the minimum number of points for coarse grained
        # ACF analysis
        nMin = 1e3
        nMin = int(nMin)
        return nMin
    
    def PathInfo(self,SystemInfo):
        SystemInfo['Path'] = {}
        SystemInfo['Path']['pwd'] = os.getcwd() + '/'
        SystemInfo['Path']['LocalRunDir'] = SystemInfo['Path']['pwd'] + 'LocalRun/'
        SystemInfo['Path']['DynamicFileDir'] = SystemInfo['Path']['pwd'] + 'DynamicFiles/'
        SystemInfo['Path']['Data'] = SystemInfo['Path']['pwd'] + '../Data/'
        
        if SystemInfo['OS'] == 'Linux':
            if SystemInfo['ComputerName'] == 'farber':
                SystemInfo['Path']['ZacrosExecuteable'] = '/home/1486/Zacros/ZacrosBinaryOutput_1.02/build/'
            else:
                print 'cluster not implemented'
        else:
            SystemInfo['Path']['ZacrosExecuteable'] = SystemInfo['Path']['pwd'] + '../Executables/'
            
        return SystemInfo