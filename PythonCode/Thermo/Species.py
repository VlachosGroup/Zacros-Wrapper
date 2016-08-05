# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 14:47:36 2016

@author: mpnun
"""

class Species(object):
    
    def __init__(self):
       
        print 'Calling Species constructor' 
        
        self.name = 'CO'
        self.elems = []
        self.Q = 1
        self.E = 0
        self.vibs = []
      
    def Q_vib(self):
        return sum(self.vibs)