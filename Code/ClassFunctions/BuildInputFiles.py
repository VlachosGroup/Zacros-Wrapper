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
        pass
    
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
            txt.write(' cd $i;qsub Zacros.ps;cd ..\n')
            txt.write('done\n')
    
    def BuildFiles(self,Cnd,BuildPath = '',SDDict = KMCut.KMCUtilities().InitializeScaleDown()):
        SystemInfo = ut.GeneralUtilities().SystemInformation()
        if BuildPath == '':
            BuildPath = SystemInfo['Path']['LocalRunDir'] + 'Build/'
        Files = ut.GeneralUtilities().GetFiles(BuildPath)
        
        #Purge Directory
        for i in Files:
            os.remove(BuildPath + i)
            
        Cnd = self.WriteSimIn(Cnd,BuildPath)
        self.WriteMechanism(Cnd,BuildPath,SDDict)
        self.WriteEnergetics(Cnd,BuildPath)
        self.WriteStateIn(Cnd,BuildPath)
        self.WriteLattice(Cnd,BuildPath)

    def WriteEnergetics(self,Cnd,BuildPath):
        nCluster = Cnd['Cluster']['nCluster']

        with open(BuildPath + 'energetics_input.dat', 'w') as txt:
            txt.write('energetics\n\n')
            for i in range(0,nCluster):
                txt.write('#'*80 + '\n\n')
                ClusterStr = Cnd['Cluster']['Input'][i]
                txt.write('cluster ' + ClusterStr['Name'] + '\n\n')
                txt.write('  sites ' + str(ClusterStr['nSites']) + '\n')
                
                if ClusterStr['neighboring'] != '':
                    txt.write('  neighboring')
                    for j in range(0,len(ClusterStr['neighboring'])):
                        txt.write(' ' + ClusterStr['neighboring'][j])
                    txt.write('\n')
                    
                txt.write('  lattice_state\n')
                for j in range(0,int(ClusterStr['nSites'])):
                    txt.write(ClusterStr['latstate'][j] + '\n')
                    
                nVariant = len(ClusterStr['variant'])              
                txt.write('\n')
                for j in range(0,nVariant):
                    txt.write('  variant ' + ClusterStr['variant'][j]['Name'] + '\n')
                    txt.write(ut.GeneralUtilities().PadStr('    site_types',25))
                    for k in range(0,len(ClusterStr['variant'][j]['site_types'])):
                        txt.write(ClusterStr['variant'][j]['site_types'][k] + ' ')
                    txt.write('\n')
                    if int(ClusterStr['variant'][j]['graph_multiplicity'])>0:
                        txt.write(ut.GeneralUtilities().PadStr('    graph_multiplicity',25) + str(ClusterStr['variant'][j]['graph_multiplicity']) + '\n')
                    txt.write(ut.GeneralUtilities().PadStr('    cluster_eng',25) + str(ClusterStr['variant'][j]['eng']) + '\n')
                    txt.write('  end_variant\n\n')
#                
                txt.write('end_cluster\n\n')
                txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_energetics')
            
    def WriteLattice(self,Cnd,BuildPath):
        with open(BuildPath + 'lattice_input.dat', 'w') as txt:
            for i in Cnd['Lattice']['Input']:
                txt.write(i + '\n')
            
    def WriteMechanism(self,Cnd,BuildPath,SDDict):
        if SDDict['Mode'] != '' and not ut.GeneralUtilities().isblank(SDDict['SDF']):
            SDBool = True
        else:
            SDBool = False
        nMech = len(Cnd['Reactions']['Input'])
        StiffCorrCounter = -1
        with open(BuildPath + 'mechanism_input.dat', 'w') as txt:
            txt.write('mechanism\n\n')
            if SDBool:
                txt.write('# Automated stiffness reconditioning employed\n')
                txt.write('# Mode: ' + SDDict['Mode'] + '\n')
                txt.write('# SDF:')
                for i in SDDict['SDF']:
                    txt.write(' ' + str(i))
                txt.write('\n\n')
            for i in range(0,nMech):
                txt.write('#'*80 + '\n\n')
                MechStr = Cnd['Reactions']['Input'][i]
                txt.write('reversible_step ' + MechStr['Name'] + '\n')
                txt.write('  sites ' + str(MechStr['nSites']) + '\n')             
                if MechStr['neighboring'] != '':
                    txt.write('  neighboring')
                    for j in range(0,len(MechStr['neighboring'])):
                        txt.write(' ' + MechStr['neighboring'][j])
                    txt.write('\n')

                if MechStr['gas_reacs_prods'] != '':
                    txt.write('  gas_reacs_prods ' + MechStr['gas_reacs_prods'][0] + ' ' + str(MechStr['gas_reacs_prods'][1]) + '\n')                    
                
                txt.write('  initial\n')
                for j in range(0,int(MechStr['nSites'])):
                    txt.write(MechStr['initial'][j] + '\n')
                    
                txt.write('  final\n')
                for j in range(0,int(MechStr['nSites'])):
                    txt.write(MechStr['final'][j] + '\n')
                    
                nVariant = len(MechStr['variant'])
                txt.write('\n')
                for j in range(nVariant):
                    txt.write('  variant ' + MechStr['variant'][j]['Name'] + '\n')
                    txt.write(ut.GeneralUtilities().PadStr('    site_types',25))
                    for k in range(0,len(MechStr['variant'][j]['site_types'])):
                        txt.write(MechStr['variant'][j]['site_types'][k] + ' ')
                    txt.write('\n')
                    pre_exp = MechStr['variant'][j]['pre_expon']
                    if SDBool:
                        StiffCorrCounter += 1
                        if Cnd['StiffnessRecondition']['APSdF'] == '':
                            APSdF = 1.
                        else:
                            APSdF = Cnd['StiffnessRecondition']['APSdF'][StiffCorrCounter]
                        rescale = SDDict['SDF'][StiffCorrCounter] / APSdF
                        pre_exp = pre_exp / rescale
                        if rescale != 1:                            
                            txt.write(ut.GeneralUtilities().PadStr('    pre_expon',25) 
                            + ut.GeneralUtilities().N2FS(pre_exp,NumType=1,digits=3))
                            txt.write('    # Pre-exponential has been rescaled by a factor of ' + 
                            ut.GeneralUtilities().N2FS(SDDict['SDF'][StiffCorrCounter],NumType=1,digits=4) 
                            + ' based on automated stiffness analysis\n')
                        else:
                            txt.write(ut.GeneralUtilities().PadStr('    pre_expon',25) 
                            + ut.GeneralUtilities().N2FS(pre_exp,NumType=1,digits=3) + '\n')
                    else:                        
                        txt.write(ut.GeneralUtilities().PadStr('    pre_expon',25) 
                        + ut.GeneralUtilities().N2FS(pre_exp,NumType=1,digits=3) + '\n')
                    txt.write(ut.GeneralUtilities().PadStr('    pe_ratio',25) + 
                        ut.GeneralUtilities().N2FS(MechStr['variant'][j]['pe_ratio'],NumType=1,digits=3) + '\n')
                    txt.write(ut.GeneralUtilities().PadStr('    activ_eng',25) + str(MechStr['variant'][j]['activ_eng']) + '\n')
                    if MechStr['variant'][j]['prox_factor'] != '':
                        txt.write(ut.GeneralUtilities().PadStr('    prox_factor',25) + str(MechStr['variant'][j]['prox_factor']) + '\n')
                    txt.write('  end_variant\n\n')
              
                txt.write('end_reversible_step\n\n')
                txt.write('#'*80 + '\n\n')
            txt.write('\n\nend_mechanism')
    
    def WriteSimIn(self,Cnd,BuildPath):
        with open(BuildPath + 'simulation_input.dat', 'w') as txt:
            SeedTxt = ''
            if Cnd['Conditions']['Seed'] == '':
                Cnd['Conditions']['Seed'] = KMCut.KMCUtilities().KMCSeed()
                SeedTxt = '      #Random seed from Python wrapper'
            
            txt.write('#KMC simulation specification\n\n')
            txt.write('random_seed' + ' '*9 + str(Cnd['Conditions']['Seed']) + SeedTxt + '\n\n')
            
            txt.write('temperature         ' + ut.GeneralUtilities().N2FS(Cnd['Conditions']['T'],NumType=3) + '\n')
            txt.write('pressure            ' + ut.GeneralUtilities().N2FS(Cnd['Conditions']['P'],NumType=3) + '\n\n') 
            txt.write('n_gas_species       ' + str(Cnd['Species']['n_gas']) + '\n')
            txt.write('gas_specs_names     ')
            for i in range(0,Cnd['Species']['n_gas']):
                txt.write(ut.GeneralUtilities().PadStr(Cnd['Species']['gas_spec'][i],14) + ' ')
                
            GasList  = ['gas_energies','gas_molec_weights','gas_molar_fracs']
            GasList2 = ['gas_eng','gas_MW','gas_molfrac']
            for j in range(0,len(GasList)):
                txt.write('\n' + ut.GeneralUtilities().PadStr(GasList[j],19) + ' ')
                for i in range(0,Cnd['Species']['n_gas']):
                    txt.write(ut.GeneralUtilities().PadStr(ut.GeneralUtilities().N2FS(Cnd['Species'][GasList2[j]][i],NumType=3),14) + ' ')

            txt.write('\n\n')
            txt.write('n_surf_species      ' + str(Cnd['Species']['n_surf']) + '\n')
            txt.write('surf_specs_names    ')
            for i in range(0,Cnd['Species']['n_surf']):
                txt.write(ut.GeneralUtilities().PadStr(Cnd['Species']['surf_spec'][i],14) + ' ')
            txt.write('\nsurf_specs_dent     ')
            for i in range(0,Cnd['Species']['n_surf']):
                txt.write(ut.GeneralUtilities().PadStr(str(Cnd['Species']['surf_dent'][i]),15))
            txt.write('\n\n')    
            
            if Cnd['Report']['hist'][0] == 'off':
                txt.write('snapshots           off\n')
            elif Cnd['Report']['hist'][0] == 'event':
                txt.write('snapshots           on '  + Cnd['Report']['hist'][0] + ' ' + str(int(Cnd['Report']['hist'][1])) + '\n')
            elif Cnd['Report']['hist'][0] == 'time':
                txt.write('snapshots           on '  + Cnd['Report']['hist'][0] + ' ' + str(np.float(Cnd['Report']['hist'][1])) + '\n')
                
            if Cnd['Report']['procstat'][0] == 'off':
                txt.write('process_statistics  off\n')
            elif Cnd['Report']['procstat'][0] == 'event':
                txt.write('process_statistics  on '  + Cnd['Report']['procstat'][0] + ' ' + str(int(Cnd['Report']['procstat'][1])) + '\n')
            elif Cnd['Report']['procstat'][0] == 'time':
                txt.write('process_statistics  on '  + Cnd['Report']['procstat'][0] + ' ' + str(np.float(Cnd['Report']['procstat'][1])) + '\n')
                
            if Cnd['Report']['specnum'][0] == 'off':
                txt.write('species_numbers     off\n')
            elif Cnd['Report']['specnum'][0] == 'event':
                txt.write('species_numbers     on '  + Cnd['Report']['specnum'][0] + ' ' + str(int(Cnd['Report']['specnum'][1])) + '\n')
            elif Cnd['Report']['specnum'][0] == 'time':
                txt.write('species_numbers     on '  + Cnd['Report']['specnum'][0] + ' ' + str(np.float(Cnd['Report']['specnum'][1])) + '\n')

            txt.write('event_report ' + Cnd['Report']['event'] + '\n\n')
            if Cnd['Conditions']['MaxStep'] == '' or re.search('inf',str(Cnd['Conditions']['MaxStep'])):
                txt.write('max_steps           infinity\n')
            else:
                txt.write('max_steps           ' + str(Cnd['Conditions']['MaxStep']) + '\n')
            
            if Cnd['Conditions']['SimTime']['Max'] == '' or re.search('inf',str(Cnd['Conditions']['SimTime']['Max'])):
                txt.write('max_time            infinity\n')
            else:
                txt.write('max_time            ' + str(Cnd['Conditions']['SimTime']['Max']) + '\n')
            txt.write('\nwall_time           ' + str(Cnd['Conditions']['WallTime']['Max']) + '\n\n')
            
            if not Cnd['Conditions']['restart']:
                txt.write('no_restart\n')
            txt.write('\nfinish\n')
            
        return Cnd

    def WriteStateIn(self,Cnd,BuildPath):
        if Cnd['StateInput']['Type'] != '':
            if Cnd['StateInput']['Type'] == 'StateInput':   #Copy from prior state_input file
                with open(BuildPath + 'state_input.dat', 'w') as txt:
                    for i in Cnd['StateInput']['Struct']:
                        txt.write(i + '\n')
            elif Cnd['StateInput']['Type'] == 'history':   #Copy from prior history_output file
                pass
            
                Lattice = Cnd['StateInput']['Struct']
                UniqSpec = np.unique(Lattice[np.not_equal(Lattice[:,2],0),1])
                nAds = len(UniqSpec)
                SpecIden = [0] * nAds
                AdsInfo = [[] for i in range(0,nAds)]
                DentInfo = [[] for i in range(0,nAds)]
                for i in range(0,nAds):
                    for j in range(0,Lattice.shape[0]):
                        if UniqSpec[i] == Lattice[j,1]:
                            AdsInfo[i].append(j+1)
                            DentInfo[i].append(Lattice[j,3])
                            SpecIden[i] = Lattice[j,2]
                
                if nAds > 0:
                    with open(BuildPath + 'state_input.dat','w') as txt:
                        txt.write('initial_state\n')
                        for i in range(0,nAds):
                            txt.write('  seed_on_sites  ' + ut.GeneralUtilities().PadStr(Cnd['Species']['surf_spec'][SpecIden[i]-1],10))
                            for j in range(0,len(DentInfo[i])):
                                for k in range(0,len(DentInfo[i])):
                                    if j + 1 == DentInfo[i][k]:
                                        txt.write(str(AdsInfo[i][k]) + '  ')
                            txt.write('\n')
                        txt.write('end_initial_state\n')
            else:
                print 'Unrecognized state_input type'
                print 'state_input not written'
                


