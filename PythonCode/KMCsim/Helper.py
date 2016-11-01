# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 15:40:49 2016

@author: robieta
"""

from itertools import (takewhile,repeat)
import numpy as np
import os
import platform
import sys

sys.path.append("..")



class Helper:
    def __init__(self):
        pass
        
    def GetFiles(self,Path):
        Files = [d for d in os.listdir(Path) if not os.path.isdir(Path + d + '/')]
        return Files
        
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
        
    def PrintDict(self,Dict,Indent = 0):
        keys = [key for key in Dict]
        keys.sort()
        for i in keys:
            if type(Dict[i]) == dict:
                if Indent == 0:
                    print ''
                print ' ' * 2 * Indent + i + ':'
                self.PrintDict(Dict[i],Indent+1)
            elif isinstance(Dict[i],np.ndarray):
                if Indent == 0:
                    print ''
                print self.PadStr(' ' * 2 * Indent + i + ':',20),
                self.PrintList(Dict[i].tolist())
            else:
                print self.PadStr(' ' * 2 * Indent + i + ':  ',20),
                self.PrintList(Dict[i])
    
    def PrintList(self,List,Nest=0):
        if len(str(List)) > 100 and len(List) > 2:
            print '[',
            if type(List[0]) == list and len(str(List[0])) > 100 and len(List[0]) > 2:
                self.PrintList(List[0],Nest+1)
            else:
                print List[0],
            print ',... ,',
            if type(List[-1]) == list and len(str(List[-1])) > 100 and len(List[-1]) > 2:
                self.PrintList(List[-1],Nest+1)
            else:
                print List[-1],
            if Nest == 0:
                print ']'
            else:
                print ']',
        else:
            print List
                
#    def SystemInformation(self):
#        SystemInfo = {}
#        SystemInfo['OS'] = platform.system()
#        if SystemInfo['OS'] == 'Linux':
#            SystemInfo['ComputerName'] = os.environ['SGE_CLUSTER_NAME']
#        else:
#            SystemInfo['ComputerName'] = platform.node()
#        SystemInfo = MS().PathInfo(SystemInfo)
#        return SystemInfo
        
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
        
    def CompareFileSize(self,Path1,Path2):
        if os.path.isfile(Path1) and os.path.isfile(Path2):        
            Size1 = os.stat(Path1).st_size
            Size2 = os.stat(Path2).st_size
            if Size1 == Size2:
                return True
            else:
                return False
        else:
            return False