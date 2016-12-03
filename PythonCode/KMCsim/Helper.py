# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 15:40:49 2016

@author: robieta
"""

from itertools import (takewhile,repeat)
import numpy as np
import os, shutil
import matplotlib.pyplot as plt
import matplotlib as mat

class Helper:

    @staticmethod
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
    
    @staticmethod    
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
    
    @staticmethod
    def PadStr(string,num):
        string = str(string)
        if len(string) < num:
            string = string + ' '*(num - len(string))
        return string
    
    @staticmethod
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
    
    @staticmethod    
    def rawbigcount(filename):
        # Taken from http://stackoverflow.com/questions/19001402/how-to-count-the-total-number-of-lines-in-a-text-file-using-python
        # Used to count the number of lines on very large files.
        with open(filename, 'rb') as f:
            bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
            nLines = sum( buf.count(b'\n') for buf in bufgen if buf )
        return nLines
        
    # Handle clearing a folder when running in parallel
    @staticmethod
    def ClearFolderContents(fldr_name):
                 
        for the_file in os.listdir(fldr_name):
            file_path = os.path.join(fldr_name, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)
    
    @staticmethod
    def PlotOptions():
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 12
                
    @staticmethod
    def PlotTrajectory(x_series, y_series, xlab = '', ylab = '', series_labels = [], fname = ''):
        
        Helper.PlotOptions()
        plt.figure()
        
        for i in range (len(y_series)):
            plt.plot(x_series[i], y_series[i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel(xlab, size=24)
        plt.ylabel(ylab, size=24)
        
        if not series_labels == []:
            plt.legend(series_labels, loc=4, prop={'size':20}, frameon=False)
        ax = plt.subplot(111)
        ax.set_position([0.2, 0.15, 0.7, 0.8])
        
        if fname == '':
            plt.show()
        else:
            plt.savefig(fname)
            plt.close()            