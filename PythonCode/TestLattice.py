# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:58:43 2016

@author: mpnun
"""

import os
from KMC_lattice import KMC_lattice
from ase.io import read
from ase.io import write 

os.system('cls')

POSCAR_fname = 'C:/Users/mpnun/Dropbox/Github/Zacros-scripts/Lattice Maker/graphene_example/graphene_template.POSCAR'

lat = KMC_lattice()
lat.workingdir = 'C:\Users\mpnun\Desktop\lat_test'
lat.mol_dat = read(POSCAR_fname)
lat.molecular_to_KMClat()
lat.Write_lattice_input()