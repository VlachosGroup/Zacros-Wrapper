# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 14:07:45 2016

@author: robieta
"""

import GeneralUtilities as ut
import os
import shutil
import subprocess
import tempfile
import time


class RunZacros:
    def __init__(self):
        pass
    
    def Run(self):
        sysinfo = ut.GeneralUtilities().SystemInformation()
        RunPath = sysinfo['Path']['LocalRunDir']
        
        try:
            Files = ut.GeneralUtilities().GetFiles(RunPath + 'Run/')
            #Purge Directory
            for i in Files:
                os.remove(RunPath + 'Run/' + i)
            
            Files = ut.GeneralUtilities().GetFiles(RunPath + 'Build/')
            for i in Files:
                shutil.move(RunPath + 'Build/' + i,RunPath + 'Run/' + i)  
        except:
            raise NameError('Failed to initialize run directory.\n' + 
                            'Check processes to see if zacros subprocess is still running')
        
        exePath = sysinfo['Path']['ZacrosExecuteable']
        
            
        if sysinfo['OS'] == 'Windows':
            with tempfile.NamedTemporaryFile() as iofile:
                p = subprocess.Popen(["cmd","/c","cd",RunPath + 'Run/','&',
                                      exePath + 'zacros_64_BinOutputExtended.exe'], 
                                     shell=False, stdout=iofile, stderr=iofile)    
            while True:
                if p.poll() is not None:
                    break
                else:
                    time.sleep(0.1)
                    
            if p.returncode != 0:
                raise NameError('Zacros exe did not run\nReturn code: ' + str(p.returncode) + 
                '\nEnsure all dependencies are installed and restart Spyder.')
        elif sysinfo['OS'] == 'Linux':
            p = subprocess.Popen("cd " + RunPath + "Run/; " + exePath + "zacros.x"
                    , stdout=subprocess.PIPE,shell=True).wait()
                            
                                              
        

        