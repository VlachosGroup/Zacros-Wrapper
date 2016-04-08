# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 15:20:36 2016

@author: robieta
"""

import GeneralUtilities as ut
import numpy as np

class ProcessOutput:
    def __init__(self):
        pass
    
    def CalcRate(self,CndIn,Stoich):
        if type(CndIn) == dict:
            CndIn = [CndIn]
            
        RateList = []
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
        print np.mean(np.array(RateList),axis=0)
        print ut.GeneralUtilities().CI(np.array(RateList),axis=0)