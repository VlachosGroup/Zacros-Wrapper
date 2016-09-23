# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 14:01:26 2016

@author: robieta
"""

import copy
import GeneralUtilities as ut
import KMCUtilities as KMCut
import numpy as np
import os
import re
import shutil

class BuildInputFiles():
    def __init__(self):
        KMCut.KMCUtilities().BuildEmptyFolders()
    
    def BuildJob(self,Cnd,SDDict = KMCut.KMCUtilities().InitializeScaleDown(),Name='',nRuns=1):
        SystemInfo = ut.GeneralUtilities().SystemInformation()
        BuildDir = SystemInfo['Path']['Data'] + 'JobBuilds/'
        if Name == '':
            Dir = ut.GeneralUtilities().GetDir(BuildDir)
            Dir = [i for i in Dir if len(i.split('_'))==2]
            for i in range(len(Dir)+1):
                if 'Job_' + str(i+1) not in Dir:
                    Name = 'Job_' + str(i+1)
                    break
        BuildDir2 = BuildDir + Name + '/'
        os.mkdir(BuildDir2)
        Cnd = KMCut.KMCUtilities().InitializeCnd(Cnd,PreserveInput=True)
        Cnd['Conditions']['Seed'] = ''
        RunList = ['0' * (len(str(nRuns))-len(str(i+1))) + str(i+1) for i in range(nRuns)]
        for i in RunList:
            os.mkdir(BuildDir2 + i + '/')
            Cnd['Conditions']['Seed'] = ''
            self.BuildFiles(Cnd,BuildPath = BuildDir2 + i + '/',SDDict=SDDict)
        SubmitScript = ['Farber_Submit.qs','Squidward_Submit.ps']
        for i in SubmitScript:
            shutil.copy(SystemInfo['Path']['pwd'] + 'SubmitFiles/' + i,BuildDir2 + i)
        self.WriteJobScript(BuildDir2,SubmitScript)
            
    def WriteJobScript(self,Path,SubmitScript):
        with open(Path + 'SubmitKMC.sh','wb') as txt:
            txt.write('#! /bin/bash\n')
            txt.write('# string to run this:\n')
            txt.write('# chmod 744 ./SubmitKMC.sh;./SubmitKMC.sh\n')
            txt.write('if [ $SGE_CLUSTER_NAME = farber ]; then\n')
            txt.write('  JobScriptName="' + SubmitScript[0] + '"\n')
            txt.write('elif [ $SGE_CLUSTER_NAME = squidward ]; then\n')
            txt.write('  JobScriptName="' + SubmitScript[1] + '"\n')
            txt.write('fi\n\n')
            txt.write('for i in $(echo ./*/);do\n')
            txt.write(' submit="${i}Zacros.ps"\n')
            txt.write(' cp ./$JobScriptName $submit\n')
            txt.write('sleep 1\n')
            txt.write(' cd $i;qsub Zacros.ps;cd ..\n')
            txt.write('done\n')