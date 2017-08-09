# Use this for converting Wei's model

import sys
import os
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
import copy
import random
from shutil import copyfile
import pickle

from ase.build import fcc111
from ase.io import read
from ase.visualize import view
from ase.io import write
from ase import Atoms


sys.path.append('/home/vlachos/mpnunez/Github/Zacros-Wrapper')
import zacros_wrapper.Lattice as lat
import zacros_wrapper as zw


'''
Create template ASE atoms objects
'''

slab_template = fcc111('Pd', size=(12, 12, 4), vacuum=15.0)


kmct = zw.kmc_traj()
kmct.path = '.'
kmct.ReadAllOutput()
kmct.KMC_lat.Read_lattice_output('lattice_output.txt')

top_lay_z = np.max(slab_template.get_positions()[:,2]) + 0.8

for snap_num in range(len(kmct.History)):

    snap = kmct.History[snap_num]
    snap_slab = copy.deepcopy(slab_template)
    
    for lat_pos in range(snap.shape[0]):
    
        site_coords = kmct.KMC_lat.cart_coords[lat_pos, :]
    
        ads = None      # vacant site
        
        if ( snap[lat_pos, 2]  == 1 ): # CO
            #ads = CO = Atoms('CO', positions = [ [site_coords[0], site_coords[1], top_lay_z], [site_coords[0], site_coords[1], top_lay_z + 1.16] ])
            ads = CO = Atoms('C', positions = [ [site_coords[0], site_coords[1], top_lay_z]])       # Use C instead of CO so that it looks different than O from the top
        elif ( snap[lat_pos, 2]  == 2 ): # O
            ads = O = Atoms('O', positions = [[site_coords[0], site_coords[1], top_lay_z]])
            
        if not ads is None:
            snap_slab.extend(ads)
            
    write(os.path.join('snapshots', 'image_' + str(snap_num + 1) + '.png'), snap_slab)