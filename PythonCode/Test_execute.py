# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 12:11:23 2016

@author: mpnun
"""

import GeneralUtilities as ut
import os
import shutil
import subprocess
import tempfile
import time

# run the executable
path = 'C:/Users/mpnun/Desktop/rescale_test/'
RunPath = path
exePath = path
OS = 'Windows'


sysinfo = ut.GeneralUtilities().SystemInformation()
RunPath2 = path

with tempfile.NamedTemporaryFile() as iofile:     
    
    print '\n--- Zacros run starting ---'   
    p = subprocess.Popen(["cmd","/c","cd",RunPath2,'&',
                          exePath + 'zacros.exe'], 
                         shell=False, stdout=iofile, stderr=iofile)    
                          
while True:
    if p.poll() is not None:
        break
    else:
        time.sleep(0.1)

if p.returncode != 0:
    raise NameError('Zacros exe did not run\nReturn code: ' + str(p.returncode) + 
    '\nEnsure all dependencies are installed and restart Spyder.')
else :
    print '--- Zacros run completed ---'                