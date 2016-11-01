# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 22:09:31 2016

@author: mpnun
"""

import numpy as np
from scipy import stats
import random
from multiprocessing import Pool

def sub_cov(mat):

    nPts = mat.shape[0]
    x_sub = np.zeros(nPts)
    y_sub = np.zeros(nPts)
    
    for j in range(nPts):
        rand_ind = random.randint(0, nPts-1)
        x_sub[j] = mat[0,rand_ind]
        y_sub[j] = mat[1,rand_ind]
    cov_mat = np.cov(x_sub,y_sub)
    return cov_mat[0,1]

class Stats:

    @staticmethod
    def mean_ci(Data,p=0.05):
        xbar = np.mean(Data)
        nPts = len(Data)
        CI = np.std(Data) * stats.t.isf(p,nPts-1)/np.sqrt(nPts)
        return [xbar, CI]
        
    @staticmethod
    def diff_ci(data1, data2, p=0.05):  # Estimate and confidence interval for the difference in sample means
        
        xbar1 = np.mean(data1)
        sd1 = np.std(data1) 
        nPts1 = len(data1)
        
        xbar2 = np.mean(data2)
        sd2 = np.std(data2) 
        nPts2 = len(data2)
        
        diff = xbar1 - xbar2
        CI = stats.t.isf(p,nPts1-1) * np.sqrt( sd1**2 / nPts1 + sd2**2 / nPts2 )
        
        return [diff, CI]
    
    @staticmethod
    def cov_ci(x, y, Nboot=100, p = 0.05, n_cores = 1):
        cov_val = Stats.cov_calc(x,y)
        
        z = np.vstack([x,y])
        if n_cores == 1:        # serial
            boot_dist = np.zeros(Nboot)
            for i in range (Nboot):           
                boot_dist[i] = sub_cov(z)
        else:               # parallel
            mats = [z for i in range(Nboot)]
            pool = Pool(processes = n_cores)
            boot_dist = pool.map(sub_cov, mats)
            pool.close()
            
        ind_high = int(round(Nboot * (1-p)) - 1)
        ind_low = int(round(Nboot * p) - 1)
        boot_dist = sorted(boot_dist)
        cov_ci = (boot_dist[ind_high] - boot_dist[ind_low]) / 2
        return [cov_val, cov_ci]
    
#    @staticmethod
#    def sub_cov(mat):
#        
#        nPts = mat.shape[0]
#        
#        x_sub = np.zeros(nPts)
#        y_sub = np.zeros(nPts)
#    
#        for j in range(nPts):
#            rand_ind = random.randint(0, nPts-1)
#            x_sub[j] = mat[0,rand_ind]
#            y_sub[j] = mat[1,rand_ind]
#        return Stats.cov_calc(x_sub, y_sub)
    
    @staticmethod
    def cov_calc(x,y):
        cov_mat = np.cov(x,y)
        return cov_mat[0,1]