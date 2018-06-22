# This file has several useful functions to be imported and used in other files

from itertools import (takewhile,repeat)
import numpy as np
import os, shutil
import matplotlib as mat
import matplotlib.pyplot as plt
import scipy.stats
import copy
from numpy.matlib import repmat


# up to 24 colors stored
colors_pool = [ 'b','yellowgreen', 'gold', 'lightcoral',
               'lightskyblue', 'darkgreen', 'orange','salmon',
               'powderblue','olivedrab', 'burlywood',  'indianred', 
               'steelblue', 'lawngreen', 'y', 'hotpink',
               'slategrey', 'palegreen', 'sandybrown', 'tomato',
               'darkviolet', 'lightgreen', 'tan','maroon']

class constant:

    ''' Gas Constants '''
    R1 = 1.9872036                 # [cal/mol-K]
    R2 = 5.189478952438986e+19     # [eV/mol-K]
    R3 = 8.6173295953479e-05       # [eV/atom-K]

    ''' Planck's Constant '''
    h1 = 6.626070040E-34           # [J-s] or [m2-kg/s]
    h2 = 1.5836687E-34             # [cal/Hz]
    h3 = 4.135667662E-15           # [eV-s]

    ''' Boltzmann's constant '''
    kb1 = 1.38064852E-23           # [J/K] or [m2-kg/s2-K]
    kb2 = 3.29982916e-24           # [cal/K]
    kb3 = 8.6173303e-5             # [eV/K]

    ''' Speed of light '''
    c1 = 186000                    # [miles/s]
    c2 = 299792458                 # [m/s]

    ''' Avogadro's number '''
    NA = 6.022140857e+23           # 1/molecules

    ''' Conversions '''
    ev_atom_2_kcal_mol = 23.06055
    amu_to_kg = 1.66053904e-27
    A2_to_m2 = 1.0e-20

    ''' Molecular Weights'''
    MW_carbon = 12.0107
    MW_hydorgen = 1.00794
    MW_oxygen = 15.9994
    MW_nitrogen = 14.0067
    

def ReadWithoutBlankLines(File,CommentLines=True):
    with open(File,'r') as Txt:
        RawTxt = Txt.readlines()
    RawTxt2 = []
    for i in RawTxt:
        if not CommentLines:
            i = i.split('#')[0]
        if len(i.split())>0:
                RawTxt2.append(i)
    return RawTxt2


def isblank(Input):
    if type(Input) != str:
        Output = False
    elif len(Input) != 1:
        Output = False
    elif Input != '':
        Output = False
    else:
        Output = True
    return Output


def PadStr(string,num):
    string = str(string)
    if len(string) < num:
        string = string + ' '*(num - len(string))
    return string


def ReturnUnique(Matrix):
    Matrix = np.array(Matrix)
    ncols = Matrix.shape[1]
    dtype = Matrix.dtype.descr * ncols
    struct = Matrix.view(dtype)
    
    uniq = np.unique(struct)
    uniq = uniq.view(Matrix.dtype).reshape(-1, ncols)
    for i in range(0,uniq.shape[0]-2):
        for j in range(i+2,uniq.shape[0]):
            if np.array_equal(uniq[i],-uniq[j]):
                uniq[j,:] = uniq[i+1,:]
                uniq[i+1,:] = -uniq[i,:]
    uniq = uniq.astype(int)  
    return uniq


def rawbigcount(filename):
    # Taken from http://stackoverflow.com/questions/19001402/how-to-count-the-total-number-of-lines-in-a-text-file-using-python
    # Used to count the number of lines on very large files.
    with open(filename, 'rb') as f:
        bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
        nLines = sum( buf.count(b'\n') for buf in bufgen if buf )
    return nLines
    
# Handle clearing a folder when running in parallel

def ClearFolderContents(fldr_name):
    
    '''
    Delete all files and folders in a directory
    '''
    
    for the_file in os.listdir(fldr_name):
        file_path = os.path.join(fldr_name, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(e)


def PlotOptions():

    '''
    Set matplotlib options
    '''

    mat.rcParams['mathtext.default'] = 'regular'
    mat.rcParams['text.latex.unicode'] = 'False'
    mat.rcParams['legend.numpoints'] = 1
    mat.rcParams['lines.linewidth'] = 2
    mat.rcParams['lines.markersize'] = 12
            

def PlotTimeSeries(x_series, y_series, xlab = 'Time (s)', ylab = '', series_labels = [], fname = '', logscale = False, spec_color_list = colors_pool):
    
    '''
    Plot multiple series against time
    '''
    
    PlotOptions()
    plt.figure()
    
    for i in range (len(y_series)):
        plt.plot(x_series[i], y_series[i], color = spec_color_list[i % len(spec_color_list)])
    
    #plt.xticks(size=20)
    #plt.yticks(size=20)
    plt.xlabel(xlab, size=24)
    plt.ylabel(ylab, size=24)
    
    if not series_labels == []:
        plt.legend(series_labels, bbox_to_anchor = (1.02,1),loc= 'upper left', prop={'size':12}, frameon=False)
    plt.tight_layout()
    
    if logscale:
        plt.yscale('log')
    
    if fname == '':
        plt.show()
    else:
        plt.savefig(fname, bbox_inches = "tight")
        plt.close()




def mean_ci(Data,p=0.05):

    '''
    Compute the mean and confidence interval (of the mean) for a set of data
    '''

    xbar = np.mean(Data)
    nPts = len(Data)
    CI = np.std(Data) * scipy.stats.t.isf(p,nPts-1)/np.sqrt(nPts)
    return [xbar, CI]
    

def diff_ci(data1, data2, p=0.05):
    
    '''
    Compute the mean of the differences and confidence interval for two sets of data
    '''
    
    xbar1 = np.mean(data1)
    sd1 = np.std(data1) 
    nPts1 = len(data1)
    
    xbar2 = np.mean(data2)
    sd2 = np.std(data2) 
    nPts2 = len(data2)
    
    diff = xbar1 - xbar2
    CI = scipy.stats.t.isf(p,nPts1-1) * np.sqrt( sd1**2 / nPts1 + sd2**2 / nPts2 )
    
    return [diff, CI]


def cov_ci(x, y, Nboot=100, p = 0.05):
    
    '''
    Compute covaraince of two data sets with bootstrapped confidence intervals
    '''
    
    B = np.vstack([np.array(x), np.array(y)])
    M = cov_mat_ci(B)
    return [M['cov_mat'][0,1], M['ci_mat'][0,1]]         

    

def cov_mat_ci(A, Nboot=100, p = 0.05):
    
    '''
    Compute full covaraince matrix with bootstrapped confidence intervals
    '''
    
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


def weighted_lin_regress(x, y, vars):

    '''
    Perform weighted linear regression on time series data
    x: vector of times
    y: vector of data
    vars: vector of variances
    
    Returns a solumn vector with the slope and intercept
    '''
    
    # Reshape data
    n_data = len(x)
    X = np.hstack([ x.reshape([n_data, 1]), np.ones([ n_data , 1 ]) ])
    Y = y.reshape([n_data, 1])
    W = np.diag(1 / vars)

    # Solve equation
    M0 = np.dot( W , X)
    M1 = np.linalg.inv( np.dot( np.transpose(X) , M0 ) )
    M2 = np.dot(M1, np.transpose(X))
    M3 = np.dot(M2, W)
    beta = np.dot(M3, Y)
    return beta