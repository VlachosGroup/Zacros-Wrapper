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
        
    def MinACFSubsample(self):
        # This function sets the minimum number of points for subsampled
        # ACF analysis
        nMin = 1e3
        nMin = int(nMin)
        return nMin
    
    def PathInfo(self,SystemInfo):
        SystemInfo['Path'] = {}
        SystemInfo['Path']['pwd'] = os.getcwd() + '/'
        SystemInfo['Path']['DynamicFileDir'] = SystemInfo['Path']['pwd']            # keep in the same directory as all the other codes
        SystemInfo['Path']['LocalRunDir'] = 'C:\Users\mpnun\OneDrive\Documents\Local_research_files\ZacrosWrapper\LocalRun/'    # where to run the stiffness reductions     
        SystemInfo['Path']['Data'] = 'C:\Users\mpnun\OneDrive\Documents\Local_research_files\ZacrosWrapper\BigJobs/'                 # where to do the long runs
        
        if SystemInfo['OS'] == 'Linux':
            if SystemInfo['ComputerName'] == 'farber':
                SystemInfo['Path']['ZacrosExecuteable'] = '../Executables/'
            else:
                print 'cluster not implemented'
        else:
            SystemInfo['Path']['ZacrosExecuteable'] = SystemInfo['Path']['pwd'] + '../Executables/'
            
        return SystemInfo