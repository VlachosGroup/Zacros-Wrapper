# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 15:40:49 2016

@author: robieta
"""

from itertools import (takewhile,repeat)
import numpy as np
import os, shutil
import sys
from mpi4py import MPI

class Helper:
    def __init__(self):
        pass
        
    def N2FS(self,num,**kwargs):     #function to convert numbers into a format that is 
                                    #suited to the fortran interpreter
        """Num Types:
        1   Floating point scientific notation (default)
        2   Integer
        3   Floating point"""
        
        NumType = 1
        digits = -1
        for key, value in kwargs.iteritems():
             if key == 'NumType':
                 NumType = value
             elif key == 'digits':
                 digits = value
        if digits > 0:
            num = np.round(num,-int(np.ceil(np.log10(num)-digits)))
        if NumType == 1:
            num = float(num)
            if num == 0:
                output = '0.0'
            elif np.log10(num) >= 0 and np.log10(num) < 1:
                output = str(num)
            else:
                output = str(num/10 ** np.floor(np.log10(num))) + 'E' +str(int(np.floor(np.log10(num))))
        elif NumType == 2:
            num = int(num)
            output = str(num)
        elif NumType == 3:
            num = np.float(num)
            output = str(num)
        return output 

    def ReadWithoutBlankLines(self,File,CommentLines=True):
        with open(File,'r') as Txt:
            RawTxt = Txt.readlines()
        RawTxt2 = []
        for i in RawTxt:
            if not CommentLines:
                i = i.split('#')[0]
            if len(i.split())>0:
                    RawTxt2.append(i)
        return RawTxt2
        
    def isblank(self,Input):
        if type(Input) != str:
            Output = False
        elif len(Input) != 1:
            Output = False
        elif Input != '':
            Output = False
        else:
            Output = True
        return Output
    
    def PadStr(self,string,num):
        string = str(string)
        if len(string) < num:
            string = string + ' '*(num - len(string))
        return string
        
    def ReturnUnique(self,Matrix):
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
        
    def rawbigcount(self,filename):
        # Taken from http://stackoverflow.com/questions/19001402/how-to-count-the-total-number-of-lines-in-a-text-file-using-python
        # Used to count the number of lines on very large files.
        with open(filename, 'rb') as f:
            bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
            nLines = sum( buf.count(b'\n') for buf in bufgen if buf )
        return nLines
        
    ''' Handle folder and file making when running in parallel'''
    @staticmethod
    def ClearFolderContents(fldr_name):
        
        COMM = MPI.COMM_WORLD
        if COMM.rank == 0:             
            for the_file in os.listdir(fldr_name):
                file_path = os.path.join(fldr_name, the_file)
                try:
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print(e) 