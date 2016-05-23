# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 15:20:36 2016

@author: robieta
"""

import GeneralUtilities as ut
import ReadOutputFiles as RO
import matplotlib.pyplot as plt
import numpy as np

class ProcessOutput:
    def __init__(self):
        pass
    
    def PlotReductionComparison(self,ID,Mode,Mode2,nSites,Stoich,PropStoich='',ColorMat='',ForceYAxis=[0,0]):
        Mean = np.zeros([len(Stoich),len(Mode),len(Mode2)])
        PropMean = np.zeros([len(Mode),len(Mode2)])
        CI = np.zeros([len(Stoich),len(Mode),len(Mode2)])
        PropCI = np.zeros([len(Mode),len(Mode2)])
        SimToWall = np.array([[0.]*len(Mode2)]*len(Mode))
        for i,M2 in enumerate(Mode2):
            for j,M in enumerate(Mode):
                CndList = RO.ReadOutputFiles().ReadJobOutput(ID + M2 + M)
                Output = self.CalcRate(CndList,nSites,Stoich,PropStoich)
                for k,S in enumerate(Stoich):
                    Mean[k,j,i] = Output['Mean'][k]
                    CI[k,j,i] = Output['CI'][k]
                if PropStoich != '':
                    PropMean[j,i] = Output['PropMean']
                    PropCI[j,i] = Output['PropCI']
                SimTime = 0
                WallTime = 0
                for Cnd in CndList:
                    SimTime += Cnd['Conditions']['SimTime']['Actual']
                    WallTime += Cnd['Conditions']['WallTime']['Actual']
                SimToWall[j,i] = SimTime/WallTime
    
        ScaleDown = (10 ** np.floor(np.log10(np.max(np.max(np.abs(
                        np.abs(Mean)+CI),axis=2),axis=1)))).tolist() 
        ScaleDown.append(10 ** np.floor(np.log10(np.max(np.max(np.abs(
                        np.abs(PropMean)+PropCI))))))
        ScaleDown = [1 for i in range(len(ScaleDown))]
        
        plt.close('all')
        nSubplot = len(Stoich) + [1 if PropStoich != '' else 0][0]
        for k in range(nSubplot):
            ax = plt.subplot(1,nSubplot,k+1)
            ax.set_xscale("log", nonposx='clip')
            for i,M2 in enumerate(Mode2):
                x=1/SimToWall[:,i]    
                if k < len(Stoich):
                    y=Mean[k,:,i]/ScaleDown[k]
                    ErrorBar=CI[k,:,i]/ScaleDown[k]
                    TitleStr = Cnd['Species']['gas_spec'][k]
                else:
                    y = PropMean[:,i]/ScaleDown[k]
                    ErrorBar=PropCI[:,i]/ScaleDown[k]
                    TitleStr = 'Propensity'

                    
                H = plt.errorbar(x,y,ErrorBar, fmt='.', capthick=2)
                if len(ColorMat) > 0:
                    plt.setp(H,'color',ColorMat[i])
            YStr = 'TOF '
            if ScaleDown[k] != 1:
                YStr += 'x 10^' + str(np.int(-np.log10(ScaleDown[k]))) + ' '
            YStr += '(1/s)'
            plt.ylabel(YStr)
            plt.xlabel('Run time / Simulation time')
            plt.title(TitleStr)
            plt.gca().set_ylim(bottom=0)
            if not np.array_equal([0,0],ForceYAxis):
                plt.ylim(np.array(ForceYAxis)/ScaleDown[k])
    
    def CalcRate(self,CndIn,nSites,Stoich,PropStoich=''):
        if type(CndIn) == dict:
            CndIn = [CndIn]
            
        RateList = []
        PropRateList = []
        WList = np.zeros((len(CndIn),len(CndIn[1]['Reactions']['Names'])))      # number of runs x number of reactions
        ind = 0
        for Cnd in CndIn:        
            #Check if run is out of burn in period
            PostBurnIn = len(['' for i in Cnd['ACF']['TauSep']['PostBurnIn'] if i != -1])>0    

            if PostBurnIn:
                nRecord = len(Cnd['Specnum']['spec'])
                BurnInTauSep = []
                for i in range(len(Cnd['ACF']['TauSep']['BurnIn'])):
                    BurnInTauSep.append(np.max(Cnd['ACF']['TauSep']['BurnIn'][i]))

                PostBurnInTauSep = np.max(Cnd['ACF']['TauSep']['PostBurnIn'])
                PostBurnInThreeTauSep = int(np.ceil(PostBurnInTauSep*3))
                BurnInTotal = int(np.sum(BurnInTauSep))
                
                FullSample = nRecord >= (BurnInTotal + PostBurnInThreeTauSep)
                if FullSample:
                    SpecnumArray = np.array(Cnd['Specnum']['spec'][BurnInTotal:])[::PostBurnInThreeTauSep,Cnd['Species']['n_surf']:]

                    RateAppend = ((SpecnumArray[1:,:]-SpecnumArray[:-1,:]) / 
                        (PostBurnInThreeTauSep * (Cnd['Specnum']['t'][1] - Cnd['Specnum']['t'][0]))/np.array(Stoich)).tolist()
                    for i in RateAppend:
                        RateList.append(i)
                        
                    if PropStoich != '':
                        PropRateTemp = np.zeros(Cnd['Binary']['propCounter'].shape[0])
                        for i,PS in enumerate(Cnd['Reactions']['Nu']):
                            if np.array_equal(PS,PropStoich):
                                PropRateTemp += Cnd['Binary']['propCounter'][:,i]
                            elif np.array_equal(PS,-np.array(PropStoich)):
                                PropRateTemp -= Cnd['Binary']['propCounter'][:,i]
                        PropRateTemp2 = PropRateTemp[BurnInTotal:][::PostBurnInThreeTauSep]
                        PropRateAppend =  ((PropRateTemp2[1:] - PropRateTemp2[:-1])/
                                (PostBurnInThreeTauSep * (Cnd['Specnum']['t'][1] - Cnd['Specnum']['t'][0])))
                        for i in PropRateAppend:
                            PropRateList.append(i)
                            
                # Store sesitiivty analysis data
                WList[ind,:] = np.array(Cnd['Binary']['W_sen_anal'][-1,:])
                ind += 1
                
            else:
                print 'Warning: Run did not exceed the burn-in period. It was not processed.'
            
        # Compute averages and confidence intervals
        Mean = np.mean(np.array(RateList),axis=0)/nSites
        PropMean = np.mean(PropRateList)/nSites
        CI = ut.GeneralUtilities().CI(np.array(RateList),axis=0)/nSites
        PropCI = ut.GeneralUtilities().CI(np.array(PropRateList),axis=0)/nSites
        
        print CndIn[1]['Reactions']['Names']
        print CndIn[1]['Reactions']['Nu']
        print CndIn[1]['Reactions']['UniqNu']       
        
        print CndIn[1]['StiffnessRecondition']['Mode']
        print CndIn[1]['StiffnessRecondition']['APSdF']
        
        # Compute terminal sensitivity values
#        print WList         # need to group these by reaction class and exclude the reactions which have been scaled down
        SenCoeff = []
        SenCoeffCI = []
        
        Output = {'Mean':Mean,'CI':CI,'PropMean':PropMean,'PropCI':PropCI,'SenCoeff':SenCoeff,'SenCoeff':SenCoeffCI}

        return Output