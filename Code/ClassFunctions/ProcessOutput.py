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
                
                print '\nSensitivity Coefficients'
                print Output['SenCoeff']
                
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
    
    def CalcRate(self,CndIn,nSites,Stoich,PropStoich=''):       # change the name of this function to reflect the fact that it only works for steady-state
        if type(CndIn) == dict:
            CndIn = [CndIn]
            
        RateList = []               # Rates computed from changes in species numbers
        PropRateList = []           # Rates computed from integral propensities, for every tau interval
        n_rxns = len(CndIn[0]['Reactions']['Names'])
        WList = np.zeros((len(CndIn),n_rxns))      # number of runs x number of reactions
        rate_obj = np.zeros((len(CndIn)))      # number of runs x number of reactions
        ind = 0
        for Cnd in CndIn:        
            
            PostBurnIn = len(['' for i in Cnd['ACF']['TauSep']['PostBurnIn'] if i != -1])>0     #Check if run is out of burn in period

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
                                rate_obj[ind] += Cnd['Binary']['prop'][-1,i]            # add the propensity at the final time to the list
                            elif np.array_equal(PS,-np.array(PropStoich)):
                                PropRateTemp -= Cnd['Binary']['propCounter'][:,i]
                                rate_obj[ind] -= Cnd['Binary']['prop'][-1,i]            # add the propensity at the final time to the list
                        PropRateTemp2 = PropRateTemp[BurnInTotal:][::PostBurnInThreeTauSep]
                        PropRateAppend =  ((PropRateTemp2[1:] - PropRateTemp2[:-1])/
                                (PostBurnInThreeTauSep * (Cnd['Specnum']['t'][1] - Cnd['Specnum']['t'][0])))
                        for i in PropRateAppend:
                            PropRateList.append(i)
                            
                # Store sesitivity analysis data
                WList[ind,:] = np.array(Cnd['Binary']['W_sen_anal'][-1,:])
                ind += 1
                
            else:
                print 'Warning: Run did not exceed the burn-in period. It was not processed.'
            
        # Compute averages and confidence intervals
        Mean = np.mean(np.array(RateList),axis=0)/nSites
        PropMean = np.mean(PropRateList)/nSites
        CI = ut.GeneralUtilities().CI(np.array(RateList),axis=0)/nSites
        PropCI = ut.GeneralUtilities().CI(np.array(PropRateList),axis=0)/nSites        
        
        # Compute sensitivity coefficients
        n_rxns_2 = n_rxns / 2
        SenCoeff = np.zeros(n_rxns_2)
        SenCoeffCI = np.zeros(n_rxns_2)
        
        print '\nTotal scaledown factor'
        print CndIn[0]['StiffnessRecondition']['APSdF']        
        
        for i in range(0,n_rxns_2):
            if CndIn[1]['StiffnessRecondition']['APSdF'][i] > 1:                # Reaction is fast and unimportant
                SenCoeff[i] = 0
                SenCoeffCI[i] = 0
            else:                                                               # Reaction is slow and may be important
                W = WList[:,2*i] + WList[:,2*i+1]                               # need to group these by reaction class
                cov_mat = np.cov(W,rate_obj)
                SenCoeff[i] = cov_mat[0,1] / np.mean(rate_obj)                  # normalize by the rate
                
                # need to add another term (1?) if the rate is the same reaction as the rate constant being perturbed                
                
                SenCoeffCI[i] = 0   # need to implement statistical bootstrapping to estimate the confidence interval
                

        Output = {'Mean':Mean,'CI':CI,'PropMean':PropMean,'PropCI':PropCI,'SenCoeff':SenCoeff,'SenCoeffCI':SenCoeffCI}

        return Output
        
        
    def CalcRateTransient(self,CndIn,nSites=1,PropStoich=''):
        if type(CndIn) == dict:
            CndIn = [CndIn]
            
        n_rxns = len(CndIn[0]['Reactions']['Names'])
        WList = np.zeros((len(CndIn),n_rxns))      # number of runs x number of reactions
        rate_obj = np.zeros((len(CndIn)))      # number of runs x number of reactions
        ind = 0
        for Cnd in CndIn:        
                
            if PropStoich != '':               
                
                for i,PS in enumerate(Cnd['Reactions']['Nu']):
                    if np.array_equal(PS,PropStoich):
                        rate_obj[ind] += Cnd['Binary']['prop'][-1,i]            # add the propensity at the final time to the list
                    elif np.array_equal(PS,-np.array(PropStoich)):
                        rate_obj[ind] -= Cnd['Binary']['prop'][-1,i]            # add the propensity at the final time to the list
            else:
                print 'Reaction rate stoichiometry not specified'
                        
            # Store sesitivity analysis data
            WList[ind,:] = np.array(Cnd['Binary']['W_sen_anal'][-1,:])
            ind += 1

            
        # Compute averages and confidence intervals
        rate_obj = rate_obj / nSites            # normalize the observable by the number of active sites
        Mean = np.mean(rate_obj)
        CI = ut.GeneralUtilities().CI(np.array(rate_obj),axis=0)          
        
        # Compute sensitivity coefficients
        n_rxns_2 = n_rxns / 2
        SenCoeff = np.zeros(n_rxns_2)
        SenCoeffCI = np.zeros(n_rxns_2)
        
        for i in range(0,n_rxns_2):
            if CndIn[0]['StiffnessRecondition']['APSdF'][i] > 1:                # Reaction is fast and unimportant
                SenCoeff[i] = 0
                SenCoeffCI[i] = 0
            else:                                                               # Reaction is slow and may be important
                W = WList[:,2*i] + WList[:,2*i+1]                               # need to group these by reaction class
                cov_mat = np.cov(W,rate_obj)
                SenCoeff[i] = cov_mat[0,1] / Mean                  # normalize by the rate
                
                # need to add another term (1?) if the rate is the same reaction as the rate constant being perturbed                
                
                SenCoeffCI[i] = 0   # need to implement statistical bootstrapping to estimate the confidence interval
                

        Output = {'Mean':Mean,'CI':CI,'SenCoeff':SenCoeff,'SenCoeffCI':SenCoeffCI}

        return Output