# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 22:09:31 2016

@author: mpnun
"""

import numpy as np
from scipy import stats
import random

class Stats:

    @staticmethod
    def mean_ci(Data,p=0.05):
        xbar = np.mean(Data)
        nPts = len(Data)
        CI = np.std(Data) * stats.t.isf(p,nPts-1)/np.sqrt(nPts)
        return [xbar, CI]    
    
    @staticmethod
    def cov_ci(x, y, Nboot=1000, p = 0.05):
        cov_val = Stats.cov_calc(x,y)
        nPts = len(x)
        
        boot_dist = np.zeros(Nboot)
        
        for i in range (Nboot):           
            x_sub = np.zeros(nPts)
            y_sub = np.zeros(nPts)
            for j in range(nPts):
                rand_ind = random.randint(0, nPts-1)
                x_sub[j] = x[rand_ind]
                y_sub[j] = y[rand_ind]
            boot_dist[i] = Stats.cov_calc(x_sub, y_sub)
        
        ind_high = int(round(Nboot * (1-p)) - 1)
        ind_low = int(round(Nboot * p) - 1)
        boot_dist = sorted(boot_dist)
        cov_ci = (boot_dist[ind_high] - boot_dist[ind_low]) / 2
        return [cov_val, cov_ci]
     
    @staticmethod
    def cov_calc(x,y):
        cov_mat = np.cov(x,y)
        return cov_mat[0,1]