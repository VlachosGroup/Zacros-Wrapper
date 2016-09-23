# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:58:43 2016

@author: mpnun
"""

# Enter data for a simple KMC lattice and then draw it.

import os
from KMC_lattice import KMC_lattice
from ase.io import read
from ase.io import write 
import numpy as np

os.system('cls')

#POSCAR_fname = 'C:/Users/mpnun/Dropbox/Github/Zacros-scripts/Lattice Maker/graphene_example/graphene_template.POSCAR'
#
#lat = KMC_lattice()
#lat.workingdir = 'C:\Users\mpnun\Desktop\lat_test'
#lat.mol_dat = read(POSCAR_fname)
#lat.molecular_to_KMClat()
#lat.Write_lattice_input()

''' Define and print lattice '''
lat = KMC_lattice()
lat.workingdir = 'C:\Users\mpnun\Desktop\lat_test'
lat.lattice_matrix = np.array([[3.0,0.0],[0.0,3.0]])
lat.site_type_names = ['top','fcc']
lat.site_type_inds = [1,1,1,1,1,1,2,2,2]
lat.frac_coords = np.array([[0.0,0.0],[0.3,0.0],[0.6,0.0],[0.0,0.3],[0.3,0.3],[0.6,0.3],[0.0,0.6],[0.3,0.6],[0.6,0.6]])
lat.neighbor_list = [[1,2],[1,4],[2,3],[2,5],[3,6],[4,5],[4,7],[5,6],[5,8],[6,9],[7,8],[8,9]]
lat.cell_list = ''          # write everything as 'self' for now
lat.Write_lattice_input()
lat.PlotLattice()