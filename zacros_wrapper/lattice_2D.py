# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:59:52 2016

@author: mpnun, Bharat Medasani
"""

import os
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from enum import Enum

class SiteType(Enum):
    TOP = 1
    T = 1
    BRIDGE = 2
    B = 2
    HOLLOW = 3
    H = 3

class Lattice2D:
    def __init__(self, matrix):
        m = np.array(matrix, dtype=np.float64).reshape((2, 2))
        lengths = np.sqrt(np.sum(m ** 2, axis=1))
        angle = np.dot(m[0], m[1]) / (lengths[0] * lengths[1])

        self._angle = np.arccos(angle) * 180. / np.pi
        self._lengths = lengths
        self._matrix = m

    def __format__(self, fmt_spec=''):
        """
        Support format printing. Supported formats are:

        1. "l" for a list format that can be easily copied and pasted, e.g.,
           ".3fl" prints something like
           "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]"
        2. "p" for lattice parameters ".1fp" prints something like
           "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}"
        3. Default will simply print a 3x3 matrix form. E.g.,
           10.000 0.000 0.000
           0.000 10.000 0.000
           0.000 0.000 10.000
        """
        m = self.matrix.tolist()
        if fmt_spec.endswith("l"):
            fmt = "[[{}, {}], [{}, {}]]"
            fmt_spec = fmt_spec[:-1]
        elif fmt_spec.endswith("p"):
            fmt = "{{{}, {}, {}}}"
            fmt_spec = fmt_spec[:-1]
            m = self.lengths_and_angle
        else:
            fmt = "{} {}\n{} {}"
        return fmt.format(*[format(c, fmt_spec) for row in m
                            for c in row])

    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix.copy())

    @property
    def matrix(self):
        """Copy of matrix representing the Lattice"""
        return np.copy(self._matrix)

    @property
    def inv_matrix(self):
        """
        Inverse of lattice matrix.
        """
        if self._inv_matrix is None:
            self._inv_matrix = np.linalg.inv(self._matrix)
        return self._inv_matrix

    def get_cartesian_coords(self, fractional_coords):
        """
        Returns the cartesian coordinates given fractional coordinates.

        Args:
            fractional_coords (3x1 array): Fractional coords.

        Returns:
            Cartesian coordinates
        """
        return np.dot(fractional_coords, self._matrix)

    def get_fractional_coords(self, cart_coords):
        """
        Returns the fractional coordinates given cartesian coordinates.

        Args:
            cart_coords (3x1 array): Cartesian coords.

        Returns:
            Fractional coordinates.
        """
        return np.dot(cart_coords, self.inv_matrix)

    @staticmethod
    def rectangular(a):
        """
        Convenience constructor for a rectangular lattice.

        Args:
            a (float): The *a* lattice parameter of the rectangular cell.

        Returns:
            Rectangular lattice of dimensions a x a.
        """
        return Lattice2D([[a, 0.0], [0.0, a]])

    @staticmethod
    def triangular(a):
        """
        Convenience constructor for a triangular lattice.

        Args:
            a (float): The *a* lattice parameter of the triangular cell.

        Returns:
            Triangular lattice of dimensions a x a.
        """
        pass

    @staticmethod
    def hexagonal(a):
        """
        Convenience constructor for a hexagonal lattice.

        Args:
            a (float): The *a* lattice parameter of the hexagonal cell.

        Returns:
            Hexagonal lattice of dimensions a x a.
        """
        pass

    @staticmethod
    def from_lengths_and_angle(ab, ang):
        """
        Create a Lattice using unit cell lengths and angle (in degrees).

        Args:
            ab (2x1 array): Lattice parameters, e.g. (4, 5).
            ang (float): Lattice angle in degrees, e.g., 90.

        Returns:
            A Lattice with the specified lattice parameters.
        """
        return Lattice2D.from_parameters(ab[0], ab[1], ang)


    @staticmethod
    def from_parameters(a, b, alpha):
        """
        Create a Lattice using unit cell lengths and angles (in degrees).

        Args:
            a (float): *a* lattice parameter.
            b (float): *b* lattice parameter.
            alpha (float): *alpha* angle in degrees.

        Returns:
            Lattice with the specified lattice parameters.
        """

        alpha_r = np.radians(alpha)

        val = np.cos(alpha_r)
        gamma_star = np.arccos(val)

        vector_a = [a, 0.0]
        vector_b = [b * np.cos(alpha_r), b * np.sin(alpha_r)]

        return Lattice2D([vector_a, vector_b])

    @classmethod
    def from_dict(cls, d):
        """
        Create a Lattice from a dictionary containing the a, b, and alpha
        parameters if fmt is None.

        Example:

            Lattice.from_dict(acell=3*[10], rprim=np.eye(3))
        """
        if "matrix" in d:
            return cls(d["matrix"])
        else:
            return cls.from_parameters(d["a"], d["b"], d["alpha"])

    @property
    def angle(self):
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        return self._angle

    @property
    def a(self):
        """
        *a* lattice parameter.
        """
        return self._lengths[0]

    @property
    def b(self):
        """
        *b* lattice parameter.
        """
        return self._lengths[1]

    @property
    def ab(self):
        """
        Lengths of the lattice vectors, i.e. (a, b)
        """
        return tuple(self._lengths)

    @property
    def alpha(self):
        """
        Angle alpha of lattice in degrees.
        """
        return self._angle

    @property
    def area(self):
        """
        Volume of the unit cell.
        """
        m = self._matrix
        return abs(np.cross(m[0], m[1]))

    @property
    def lengths_and_angle(self):
        """
        Returns (lattice lengths, lattice angles).
        """
        return tuple(self._lengths), self._angle

    def __repr__(self):
        outs = ["Lattice", "    ab : " + " ".join(map(repr, self._lengths)),
                " angles : " + " ".join(map(repr, self._angles)),
                " area : " + repr(self.area),
                "      A : " + " ".join(map(repr, self._matrix[0])),
                "      B : " + " ".join(map(repr, self._matrix[1]))]
        return "\n".join(outs)

    def __eq__(self, other):
        """
        A lattice is considered to be equal to another if the internal matrix
        representation satisfies np.allclose(matrix1, matrix2) to be True.
        """
        if other is None:
            return False
        # shortcut the np.allclose if the memory addresses are the same
        # (very common in Structure.from_sites)
        return self is other or np.allclose(self.matrix, other.matrix)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return "\n".join([" ".join(["%.6f" % i for i in row])
                          for row in self._matrix])

    def as_dict(self, verbosity=0):
        """""
        Json-serialization dict representation of the Lattice.

        Args:
            verbosity (int): Verbosity level. Default of 0 only includes the
                matrix representation. Set to 1 for more details.
        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "matrix": self._matrix.tolist()}
        if verbosity > 0:
            d.update({
                "a": float(self.a),
                "b": float(self.b),
                "alpha": float(self.alpha),
                "area": float(self.area)
            })

        return d

    def get_distance_and_image(self, frac_coords1, frac_coords2, jimage=None):
        """
        Gets distance between two frac_coords assuming periodic boundary
        conditions. If the index jimage is not specified it selects the j
        image nearest to the i atom and returns the distance and jimage
        indices in terms of lattice vector translations. If the index jimage
        is specified it returns the distance between the frac_coords1 and
        the specified jimage of frac_coords2, and the given jimage is also
        returned.

        Args:
            fcoords1 (2x1 array): Reference fcoords to get distance from.
            fcoords2 (2x1 array): fcoords to get distance from.
            jimage (2x1 array): Specific periodic image in terms of
                lattice translations, e.g., [1,0] implies to take periodic
                image that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            (distance, jimage): distance and periodic lattice translations
            of the other site for which the distance applies. This means that
            the distance between frac_coords1 and (jimage + frac_coords2) is
            equal to distance.
        """
        if jimage is None:
            v, d2 = pbc_shortest_vectors(self, frac_coords1, frac_coords2,
                                         return_d2=True)
            fc = self.get_fractional_coords(v[0][0]) + frac_coords1 - \
                 frac_coords2
            fc = np.array(np.round(fc), dtype=np.int)
            return np.sqrt(d2[0, 0]), fc

        mapped_vec = self.get_cartesian_coords(jimage + frac_coords2
                                               - frac_coords1)
        return np.linalg.norm(mapped_vec), jimage


class Site2D:
    pass


class ZachrosLattice:
    '''Handles the KMC lattice'''
    fname_in = 'lattice_input.dat'
    fname_out = 'lattice_output.txt'
    
    def __init__(self, lattice, size, sites, connectivity):
        ''' Initialize class variables '''
        
        #self.text_only = True  # If true, only store the text from the input file, if false, store all the complex information
        #self.lattice_in_txt = None
        
        self._matrix = matrix         # each row is a lattice vector
        self.repeat = [1,1]
        self.site_type_names = []
        self.site_type_inds = []
        self.frac_coords = []
        self.cart_coords = []
        self.neighbor_list = []
        self.cell_list = [] # self, north, northeast, east, or southeast

        
    def read_input(self, fldr):
        '''Read lattice_input.dat
        
        :param fldr: Folder directory from which to read lattice_input.dat
        '''
        self.lattice_in_txt = []
        with open(os.path.join(fldr, self.fname_in),'r') as Txt:
            RawTxt = Txt.readlines()   
        for i in RawTxt:
            self.lattice_in_txt.append(i.split('\n')[0])

    def write_input(self, fldr):
        '''Write lattice_input.dat
        
        :param fldr: Folder directory from which to read lattice_input.dat
        '''
        if self.text_only:
            with open(os.path.join(fldr, self.fname_in), 'w') as txt:
                for i in self.lattice_in_txt:
                    txt.write(i + '\n')
        else:
            self.Write_lattice_input(fldr)
            
    def set_frac_coords(self, fc):
        '''Set the fractional and Cartesian coordinates
        
        :param fc: n x 3 array (or list of lists) of fractional
            coordinates to set for each atom
        '''
        self.frac_coords = np.array(fc)
        self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)

    def set_cart_coords(self, cc):
        '''Set the Cartesian and fractional coordinates
        
        :param fc: n x 3 array (or list of lists) of Cartesian
            coordinates to set for each atom
        '''
        self.cart_coords = np.array(cc)
        self.frac_coords = np.dot(self.cart_coords,
                                  np.linalg.inv(self.lattice_matrix))
        
    def coord_shift(self, a_ind, b_ind):
        '''
        Give the coordinates of the periodic image of B which is
        closest to A
        
        :param a_ind: Index of atom A
        
        :param b_ind: Index of atom B
        '''
        a_coords = self.cart_coords[ a_ind , : ]

        image_coords = [ self.cart_coords[ b_ind , : ] for i in range(9) ]
        dist_list = [None for i in range(9)]
        ind = 0
        for we in [-1, 0, 1]:
            for sn in [-1, 0, 1]:
                image_coords[ind] += np.dot(np.array([we, sn]),
                                            self.lattice_matrix)
                dist_list[ind] = np.linalg.norm(image_coords[ind] - a_coords)
                ind += 1

        min_d = dist_list[0]
        min_d_ind = 0
        for i in range(1,9):
            if dist_list[i] < min_d:
                min_d = dist_list[i]
                min_d_ind = i
        
        return image_coords[min_d_ind]
    
    def plot(self, cutoff=3.0, plot_neighbs=False,
             type_symbols=['o','s','^','v', '<', '>', '8', 'd', 'D', 'H',
                           'h', '*', 'p', '+', ',', '.', '1', '2', '3', '4',
                           '_', 'x', '|', 0, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8],
             ms=7):
        '''
        :param cutoff: Maximum distance to draw connections between
            nearest neighbor sites. This prevents drawing line segments
            between sites which are neighbors only though their periodic
            images.
        
        :param plot_neighbs: Flag to plot line segments between lattice
            sites which are first nearest neighbors to each other
        
        :param type_symbols: List of symbols for each lattice site type.
        
        :param ms: Marker size
        
        :returns: pyplot object with the lattice graphed on it
        '''
        if self.text_only:
            raise NameError('Lattice must be built before plotting. '
                            'Use ReadAllOutput(build_lattice=True)')
        
        if self.cart_coords == []:
            self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
        
        border = np.dot(np.array(
            [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]),
            self.lattice_matrix)
        
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 12
        
        plt.figure()
        
        plt.plot(border[:,0], border[:,1], '--k', linewidth = 2) # cell border

        if plot_neighbs:
            ind = 0
            for pair in self.neighbor_list: # neighbors
                p1 = self.cart_coords[pair[0],:]
                p2 = self.cart_coords[pair[1],:]
                if self.cell_list[ind] == 'self':
                    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], '-k', linewidth=1)
                ind += 1

        for site_type in range(1, np.max(np.array(self.site_type_inds))+1):

            is_of_type = []
            
            for site_ind in range(len(self.site_type_inds)):
                if self.site_type_inds[site_ind] == site_type:
                    is_of_type.append(site_ind)
            
            plt.plot(self.cart_coords[is_of_type,0],
                     self.cart_coords[is_of_type,1], linestyle='None',
                     marker=type_symbols[(site_type-1) % len(type_symbols)],
                     color=[0.9, 0.9, 0.9], markersize=ms,
                     markeredgewidth=0.0)          # sites
        
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
        plt.xlabel('x-coord (ang)', size=24)
        plt.ylabel('y-coord (ang)', size=24)
        plt.axis('equal')
        plt.tight_layout()
        
        return plt
        
    def read_lattice_output(self, fldr):
        '''Read lattice_output.txt'''

        with open( os.path.join(fldr, self.fname_out), 'r') as txt:
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
                if int(site_2) > 0:    # Zeros are placeholders in the output file
                    self.neighbor_list.append([site_ind+1, int(site_2)])
        
        # Convert to fractional coordinates
        self.frac_coords = np.dot(self.cart_coords,
                                  np.linalg.inv(self.lattice_matrix))
    
        self.text_only = False

    def write_lattice_input(self, fldr):
        '''Write lattice_input.dat'''
    
        with open(os.path.join(fldr, self.fname_in), 'w') as txt:
            txt.write('# Lattice specification file: generated by '
                      'ZacrosWrapper' + '\n\n')
            txt.write('lattice periodic_cell\n\n');
            txt.write('cell_vectors       # in row format (Angstroms)\n')
            txt.write('\t {0:.3f} \t'.format(self.lattice_matrix[0,0]) )
            txt.write('{0:.3f} \n'.format(self.lattice_matrix[0,1]) )
            txt.write('\t {0:.3f} \t'.format(self.lattice_matrix[1,0]) )
            txt.write('{0:.3f} \n\n'.format(self.lattice_matrix[1,1]) )

            txt.write('repeat_cell\t {} \t {} \n\n'.format(self.repeat[0],
                                                           self.repeat[1]))

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
            txt.write('site_coordinates \t # fractional coordinates (x,y) in '
                      'row format\n')
            for i in range(0, len(self.site_type_inds)):
                txt.write('\t {0:.3f} \t'.format(self.frac_coords[i,0]))
                txt.write('{0:.3f} \n'.format(self.frac_coords[i,1]))
            txt.write('\n')
            
            # Site neighboring structure
            txt.write('neighboring_structure \t # site-neighsite cell\n')
            ind = 0
            for pair in self.neighbor_list:
                txt.write('\t {}-{} {} \n'.format(pair[0]+1,pair[1]+1,
                                                  self.cell_list[ind]))
                ind += 1
            txt.write('end_neighboring_structure\n\n')
            
            txt.write('end_lattice\n')
        txt.close()
        
    def build_neighbor_list(self, cut = 3.0, cut_mat = []):
        # appending lists causes this method to be slow
        '''
        Builds the neighbor list based on distances between sites
        
        :param cut: Maximum distance between nearest neighbor sites or
            their periodic images.
        
        :param cut_mat: List of nearest neighbor cutoff distances
            between different site types. If empty, it will assume the
            same distance for all site types.
        '''
        # Make sure these lists are empty before we start appending
        self.neighbor_list = []
        self.cell_list = [] # self, north, northeast, east, or southeast
    
        if self.cart_coords == []:
            self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
    
        n_site_types = len(self.site_type_names)
        n_sites = len(self.site_type_inds)
        self.cart_coords = np.dot(self.frac_coords, self.lattice_matrix)
    
        if cut_mat == []: # Use a matrix to have different cutoff distances for different site types
            cut_mat = cut * np.ones([n_site_types, n_site_types])
            
        # Loop through all sites
        for site_1 in range(n_sites):
            for site_2 in range(n_sites):
            
                c1 = self.cart_coords[site_1,:]
                c2 = self.cart_coords[site_2,:]
            
                c_2_north = c2 + self.lattice_matrix[1,:]
                c_2_northeast = c2 + self.lattice_matrix[0,:] + \
                                self.lattice_matrix[1,:]
                c_2_east = c2 + self.lattice_matrix[0,:]
                c_2_southeast = c2 + self.lattice_matrix[0,:] - \
                                self.lattice_matrix[1,:]
            
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
        
