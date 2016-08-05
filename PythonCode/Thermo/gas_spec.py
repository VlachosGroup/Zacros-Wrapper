# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 14:54:51 2016

@author: mpnun
"""

from Species import Species

class gas_spec(Species):
    
    def __init__(self):
        super(gas_spec, self).__init__()
        print 'Calling surf_spec constructor'
      
    def method(self):
        print 'running method'