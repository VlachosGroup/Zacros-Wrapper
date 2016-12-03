# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 22:09:31 2016

@author: mpnun
"""

import numpy as np
from scipy import stats

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
    def cov_ci(x, y, Nboot=100, p = 0.05):
        
        cov_val = Stats.cov_calc(x,y)

        n_points = len(x)
        boot_dist = np.zeros(Nboot)
        for i in range (Nboot):
            subpop_inds = np.random.randint(n_points, size=n_points)
            x_sub = x[subpop_inds]
            y_sub = y[subpop_inds]
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
        
    @staticmethod
    def cov_mat_ci(A, Nboot=100, p = 0.05):
        
        x = A.shape
        n_vars = x[0]
        n_obs = x[1]
        
        pop = np.zeros([n_vars, n_vars, Nboot])
        
        # Compute distribution of covariance estimates
        for i in range(Nboot):
            subpop_inds = np.random.randint(n_obs, size=n_obs)
            pop[:,:,i] = np.cov(A[:,subpop_inds])
        
        # Sort covariance estimates
        for var1 in range(n_vars):
            for var2 in range(n_vars):
                pop[var1,var2,:] = sorted(pop[var1,var2,:])
        
        # Compute half-lengths of the confidence intervals
        ind_high = int(round(Nboot * (1-p)) - 1)
        ind_low = int(round(Nboot * p) - 1)
        ci_mat = ( pop[:,:,ind_high] - pop[:,:,ind_low]) / 2.0        
        
        return {'cov_mat': np.cov(A), 'ci_mat': ci_mat}