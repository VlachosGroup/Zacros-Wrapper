# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:59:52 2016

@author: mpnun
"""

import os
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt

class Lattice:
    
    '''
    Class which handles the KMC lattice
    '''
    
    fname_in = 'lattice_input.dat'
    fname_out = 'lattice_output.txt'
    
    def __init__(self):
        
        '''
        Initialize class variables
        '''
        
        self.text_only = True                    # If true, only store the text from the input file, if false, store all the complex information
        self.lattice_in_txt = None
        
        self.lattice_matrix = None         # each row is a lattice vector
        self.repeat = [1,1]
        self.site_type_names = []
        self.site_type_inds = []
        self.frac_coords = []
        self.cart_coords = []
        self.neighbor_list = []
        self.cell_list = []     # self, north, northeast, east, or southeast

        
    def ReadIn(self, fldr):
    
        '''
        Read lattice_input.dat
        '''
    
        self.lattice_in_txt = []
        with open(os.path.join(fldr, self.fname_in),'r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            self.lattice_in_txt.append(i.split('\n')[0])

    
    def WriteIn(self, fldr):
    
        '''
        Write lattice_input.dat
        '''
        
        if self.text_only:
        
            with open(os.path.join(fldr, self.fname_in), 'w') as txt:
                for i in self.lattice_in_txt:
                    txt.write(i + '\n')
                    
        else:
        
            self.Write_lattice_input(fldr)
            
    
    def set_frac_coords(self, fc):
        
        '''
        Set the fractional and Cartesian coordinates
        '''
        
        self.frac_coords = np.array(fc)
        self.cart_coords = self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
        
        
    def set_cart_coords(self, cc):
    
        '''
        Set the Cartesian and fractional coordinates
        '''
    
        self.cart_coords = np.array(cc)
        self.frac_coords = np.dot(self.cart_coords, np.linalg.inv(self.lattice_matrix))
        

    def coord_shift(self, a_ind, b_ind):

        '''
        Give the coordinates of the periodic image of B which is closest to A
        '''

        a_coords = self.cart_coords[ a_ind , : ]

        image_coords = [ self.cart_coords[ b_ind , : ] for i in range(9) ]
        dist_list = [None for i in range(9)]
        ind = 0
        for we in [-1, 0, 1]:
            for sn in [-1, 0, 1]:
                image_coords[ind] = image_coords[ind] + np.dot(np.array([we, sn]), self.lattice_matrix)
                dist_list[ind] = np.linalg.norm( image_coords[ind] - a_coords )
                ind += 1

        min_d = dist_list[0]
        min_d_ind = 0
        for i in range(1,9):
            if dist_list[i] < min_d:
                min_d = dist_list[i]
                min_d_ind = i
        
        return image_coords[min_d_ind]
    
    
    def PlotLattice(self, cutoff = 3.0, plot_neighbs = False, type_symbols = ['o','s','^','v', '<', '>', '8', 
        'd', 'D', 'H', 'h', '*', 'p', '+', ',', '.', '1', '2', '3', '4', '_', 'x', '|', 0, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8], ms = 4):
        
        '''
        Return a pyplot object with the lattice graphed on it
        '''
        
        if self.cart_coords == []:
            self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
        
        border = np.dot(np.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]), self.lattice_matrix)             
        
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 12
        
        plt.figure()
        
        plt.plot(border[:,0], border[:,1], '--k', linewidth = 2)                  # cell border 

        if plot_neighbs:
            ind = 0
            for pair in self.neighbor_list:                                             # neighbors
                p1 = self.cart_coords[pair[0],:]
                p2 = self.cart_coords[pair[1],:]
                if self.cell_list[ind] == 'self':
                    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], '-k', linewidth = 1)
                ind += 1

        for site_type in range(1, np.max(np.array(self.site_type_inds))+1 ):

            is_of_type = []
            
            for site_ind in range(len(self.site_type_inds)):
                if self.site_type_inds[site_ind] == site_type:
                    is_of_type.append(site_ind)
            
            plt.plot(self.cart_coords[is_of_type,0], self.cart_coords[is_of_type,1], linestyle='None', marker = type_symbols[(site_type-1) % len(type_symbols)], color = [0.8, 0.8, 0.8], markersize = ms)          # sites  
        
        # Choose range to plot
        xmin = np.min(border[:,0])
        xmax = np.max(border[:,0])
        delx = xmax - xmin
        ymin = np.min(border[:,1])
        ymax = np.max(border[:,1])
        dely = ymax - ymin
        mag = 0.1        
        
        plt.xlim([xmin - mag * delx, xmax + mag * delx])
        plt.ylim([ymin - mag * dely, ymax + mag * dely])
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('x-coord (ang)',size=24)
        plt.ylabel('y-coord (ang)',size=24)
        plt.axis('equal')
        plt.tight_layout()
        
        return plt
        
    
    def Read_lattice_output(self, fldr):
    
        '''
        Read lattice_output.txt
        '''
    
        with open( os.path.join(fldr, self.fname_out) ,'r') as txt:
            RawTxt = txt.readlines()           
        n_sites = len(RawTxt) - 2
        self.cart_coords = np.zeros([n_sites,2])        
        
        # Fill in lattice vectors
        self.lattice_matrix = np.zeros([2,2])
        line1 = RawTxt[0].split()
        self.lattice_matrix[0,0] = float(line1[1])
        self.lattice_matrix[0,1] = float(line1[2])
        line2 = RawTxt[1].split()
        self.lattice_matrix[1,0] = float(line2[1])
        self.lattice_matrix[1,1] = float(line2[2])
        
        # Fill in site coordinates and neighbors
        self.neighbor_list = []
        self.site_type_inds = []
        for site_ind in range(n_sites):
            line = RawTxt[site_ind+2].split()
            self.cart_coords[site_ind,:] = [line[1], line[2]]
            self.site_type_inds.append(int(line[3]))
            neighbs = line[5::]
            for site_2 in neighbs:
                if int(site_2) > 0:         # Zeros are placeholders in the output file
                    self.neighbor_list.append([site_ind+1, int(site_2)])
        
        # Convert to fractional coordinates
        self.frac_coords = np.dot(self.cart_coords, np.linalg.inv(self.lattice_matrix))
    
        self.text_only = False

    
    
    def Write_lattice_input(self, fldr):

        '''
        Write lattice_input.dat
        '''
    
        with open(os.path.join(fldr, self.fname_in), 'w') as txt:
            txt.write('# Lattice specification file: generated by ZacrosWrapper' + '\n\n')
            txt.write('lattice periodic_cell\n\n');
            txt.write('cell_vectors       # in row format (Angstroms)\n')
            txt.write('\t {0:.3f} \t'.format(self.lattice_matrix[0,0]) )
            txt.write('{0:.3f} \n'.format(self.lattice_matrix[0,1]) )
            txt.write('\t {0:.3f} \t'.format(self.lattice_matrix[1,0]) )
            txt.write('{0:.3f} \n\n'.format(self.lattice_matrix[1,1]) )

            txt.write('repeat_cell\t {} \t {} \n\n'.format(self.repeat[0],self.repeat[1]))

            txt.write('n_cell_sites \t {} \n'.format(len(self.site_type_inds)))
            txt.write('n_site_types \t {} \n'.format(len(self.site_type_names)))
            
            txt.write('site_type_names \t')
            for site_type in self.site_type_names:
                txt.write(site_type + '\t')
            txt.write('\n')
            
            txt.write('site_types \t')
            for site_type_ind in self.site_type_inds:
                txt.write('{}  '.format(site_type_ind))
            txt.write('\n\n')
            
            # Site coordinates
            txt.write('site_coordinates \t # fractional coordinates (x,y) in row format\n')
            for i in range(0, len(self.site_type_inds)):
                txt.write('\t {0:.3f} \t'.format(self.frac_coords[i,0]))
                txt.write('{0:.3f} \n'.format(self.frac_coords[i,1]))
            txt.write('\n')
            
            # Site neighboring structure
            txt.write('neighboring_structure \t # site-neighsite cell\n');
            ind = 0
            for pair in self.neighbor_list:
                txt.write('\t {}-{} {} \n'.format(pair[0]+1,pair[1]+1,self.cell_list[ind]))
                ind += 1
            txt.write('end_neighboring_structure\n\n')
            
            txt.write('end_lattice\n')
        txt.close()
        
        
    def Build_neighbor_list(self, cut = 3.0, cut_mat = []):     # appending lists causes this method to be slow
    
        '''
        Builds the neighbor list based on distances between sites
        '''
    
        # Make sure these lists are empty before we start appending
        self.neighbor_list = []
        self.cell_list = []     # self, north, northeast, east, or southeast
    
        if self.cart_coords == []:
            self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
    
        n_site_types = len(self.site_type_names)
        n_sites = len(self.site_type_inds)
        self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
    
        if cut_mat == []:
            cut_mat = cut * np.ones([n_site_types, n_site_types])        # Use a matrix to have different cutoff distances for different site types
            
        # Loop through all sites
        for site_1 in range(n_sites):
            for site_2 in range(n_sites):
            
                c1 = self.cart_coords[site_1,:]
                c2 = self.cart_coords[site_2,:]
            
                c_2_north = c2 + self.lattice_matrix[1,:]
                c_2_northeast = c2 + self.lattice_matrix[0,:] + self.lattice_matrix[1,:]
                c_2_east = c2 + self.lattice_matrix[0,:]
                c_2_southeast = c2 + self.lattice_matrix[0,:] - self.lattice_matrix[1,:]
            
                if site_1 < site_2:        # check self
                    if np.linalg.norm( c1 - c2 ) < cut:
                        self.neighbor_list.append([site_1, site_2])
                        self.cell_list.append('self')
                        
                if np.linalg.norm( c1 - c_2_north ) < cut:
                    self.neighbor_list.append([site_1, site_2])
                    self.cell_list.append('north')
                        
                if np.linalg.norm( c1 - c_2_northeast ) < cut:
                    self.neighbor_list.append([site_1, site_2])
                    self.cell_list.append('northeast')
                    
                if np.linalg.norm( c1 - c_2_east ) < cut:
                    self.neighbor_list.append([site_1, site_2])
                    self.cell_list.append('east')
                    
                if np.linalg.norm( c1 - c_2_southeast ) < cut:
                    self.neighbor_list.append([site_1, site_2])
                    self.cell_list.append('southeast')
        