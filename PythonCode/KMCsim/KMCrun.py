# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 14:54:26 2016

@author: robieta
"""

from KMCrun_data import KMCrun_data
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat

# For executable
import GeneralUtilities as ut
import os
import signal
import shutil
import subprocess
import tempfile
import time
import matplotlib.animation as animation

class KMCrun:
    
    def __init__(self):
        
        self.data = KMCrun_data()
        self.exe_file = ''
        self.op_system = 'Windows'
        self.anim = []          # animation object used for lattice movie
        
    def Run_sim(self):
        
        os.chdir(self.data.Path)
        
        print '--- Zacros run starting ---'
        
        if self.op_system == 'Windows':
            subprocess.call([self.exe_file])
        elif self.op_system == 'Linux':
            print 'Linux execution to be implemented'
        else:
            print 'Unknown operating system'
        
        print '--- Zacros run completed ---'

    def PlotOptions(self):
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 12

    def PlotSurfSpecVsTime(self):
        self.PlotOptions()
        plt.figure()
        
        for i in range (len(self.data.Species['surf_spec'])):
            plt.plot(self.data.Specnum['t'], self.data.Specnum['spec'][:,i])    
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('spec. pop.',size=24)
#        plt.ylabel('coverage',size=24)
        plt.legend(self.data.Species['surf_spec'],loc=4,prop={'size':20},frameon=False)        
        plt.show()
    
    def PlotGasSpecVsTime(self):
        self.PlotOptions()
        plt.figure()          
          
        for i in range (len(self.data.Species['gas_spec'])):
            ind = i + len(self.data.Species['surf_spec'])
            plt.plot(self.data.Specnum['t'], self.data.Specnum['spec'][:,ind])    
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('spec. pop.',size=24)
        plt.legend(self.data.Species['gas_spec'],loc=2,prop={'size':20},frameon=False)        
        plt.show()     
    
    def PlotPropsVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.data.Reactions['names'])):
            if np.max(np.abs(self.data.Binary['prop'][:,i])) > 0:
                plt.plot(self.data.Specnum['t'], self.data.Binary['prop'][:,i]) 
                labels.append(self.data.Reactions['names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('props',size=24)
        plt.legend(labels,loc=2,prop={'size':20},frameon=False)        
        plt.show()
        
    def PlotIntPropsVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.data.Reactions['names'])):
            if np.max(np.abs(self.data.Binary['propCounter'][:,i])) > 0:
                plt.plot(self.data.Specnum['t'], self.data.Binary['propCounter'][:,i]) 
                labels.append(self.data.Reactions['names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('integral props',size=24)
        plt.legend(labels,loc=2,prop={'size':20},frameon=False)        
        plt.show()
    
    def PlotWVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.data.Reactions['names'])):
            if np.max(np.abs(self.data.Binary['W_sen_anal'][:,i])) > 0:
                plt.plot(self.data.Specnum['t'], self.data.Binary['W_sen_anal'][:,i]) 
#                plt.plot(self.data.Specnum['t'], self.data.Procstat['events'][:,i] - self.data.Binary['propCounter'][:,i]) 
                labels.append(self.data.Reactions['names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('W',size=24)
        plt.legend(labels, loc=2,prop={'size':20},frameon=False)        
        plt.show()       
    
    def PlotElemStepFreqs(self):
        self.PlotOptions
        plt.figure()        
        
        width = 0.2
        ind = 0
        yvals = []
        ylabels = []

        for i in range (self.data.Reactions['nrxns']):
            if self.data.Procstat['events'][-1,2*i] + self.data.Procstat['events'][-1,2*i+1] > 0:
                net_freq = abs(self.data.Procstat['events'][-1,2*i] - self.data.Procstat['events'][-1,2*i+1])
                if self.data.Procstat['events'][-1,2*i] > 0:              
                    plt.barh(ind-0.4, self.data.Procstat['events'][-1,2*i], width, color='r')
                if self.data.Procstat['events'][-1,2*i+1] > 0:
                    plt.barh(ind-0.6, self.data.Procstat['events'][-1,2*i+1], width, color='b')
                if net_freq > 0:
                    plt.barh(ind-0.8, net_freq, width, color='g')
                ylabels.append(self.data.Reactions['names'][i])
                yvals.append(ind-0.6)
                ind = ind - 1

        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('frequency',size=24)
        plt.xscale('log')
        plt.yticks(yvals, ylabels)
        plt.legend(['fwd','bwd','net'],loc=4,prop={'size':20},frameon=False)
        ax = plt.subplot(111)        
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        plt.show()
    
    def LatticeMovie(self):       # Need to complete this function by plotting adsorbates from the history file data

        self.data.KMC_lat.Read_lattice_output(self.data.Path + 'lattice_output.txt')
        cart_coords = self.data.KMC_lat.cart_coords
        lat = self.data.KMC_lat.lattice_matrix
#        self.data.KMC_lat.PlotLattice()
        border = np.dot(np.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]), lat)
        
        self.PlotOptions

        fig = plt.figure()
        ax = plt.axes(xlim=[np.min(border[:,0]), np.max(border[:,0])], ylim=[np.min(border[:,1]), np.max(border[:,1])])
        plt.axis('equal')        
        
        # Plot cell border 
        plt.plot(border[:,0], border[:,1], '--k', linewidth = 4)
        
        # Plot neighbors
        cutoff = 3.0
        for pair in self.data.KMC_lat.neighbor_list:
                p1 = np.array([cart_coords[pair[0]-1,0], cart_coords[pair[0]-1,1]])
                p2 = np.array([cart_coords[pair[1]-1,0], cart_coords[pair[1]-1,1]])
                if np.linalg.norm(p2 - p1) < cutoff:
                    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], '-k', linewidth = 1)
        
        # Plot sites            
        plt.plot(cart_coords[:,0], cart_coords[:,1], 'o', markersize = 15)
        
        line, = ax.plot([], [], 'rs')
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('x-coord (ang)',size=24)
        plt.ylabel('y-coord (ang)',size=24)        
        
        
        # initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,
        
        # animation function.  This is called sequentially
        def animate(i):
            
            x_list = []
            y_list = []
            
            snap = self.data.History[i]            
            for site_ind in range(self.data.n_sites):
                if snap[site_ind,2] > 0:
                    x_list.append(cart_coords[site_ind,0])
                    y_list.append(cart_coords[site_ind,1])
        
            x = np.array(x_list)
            y = np.array(y_list)
            
            line.set_data(x, y)
            return line,
        
        # call the animator.  blit=True means only re-draw the parts that have changed.
        self.anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=self.data.n_snapshots, interval=20, blit=True)
        
        plt.show()
    
    def ComputeTOF(self,Product):                       # return TOF and TOF error
        
        # Find the index of the product species
        product_ind = -1       
        for i in enumerate(self.data.Species['gas_spec']):
            if i[1] == Product:
                product_ind = i[0]
        
        # Make sure the index has been found
        if product_ind == -1:
            print 'Product species not found'
        else:
            product_ind = product_ind + self.data.Species['n_surf']         # Adjust index to account for surface species   
        
        
        nRxns = len(self.data.Reactions['Nu'])        
        TOF_contributions = [0 for i in range(nRxns)]              # number of product molecules produced in each reaction        
        for i, elem_stoich in enumerate(self.data.Reactions['Nu']):
            TOF_stoich = elem_stoich[product_ind]
            r = self.data.Binary['propCounter'][-1,i] / self.data.Specnum['t'][-1]      # ergodic average
#            r = self.data.Binary['prop'][-1,i]                                           # non-ergodic average
            TOF_contributions[i] = TOF_stoich * r         
               
        TOF = np.sum(TOF_contributions)
        TOF_fracs = TOF_contributions / TOF             # will need this for sensitivity analysis
#        return TOF
        return {'TOF': TOF, 'TOF_fracs': TOF_fracs}    
    
    def AdjustPreExponentials(self, delta_sdf):
        
        rxn_ind = 0
        for rxn_type in self.data.Reactions['Input']:
            for variant in rxn_type['variant']:
                variant['pre_expon'] = variant['pre_expon'] * delta_sdf[rxn_ind]
                self.data.scaledown_factors[rxn_ind] = self.data.scaledown_factors[rxn_ind] * delta_sdf[rxn_ind]
                rxn_ind += 1    
    
    def CheckSteadyState(self, Product, frac_sample = 0.2, d_cut = 0.12, show_graph = False):
        
        # Find the index of the product species
        product_ind = -1       
        for spec in enumerate(self.data.Species['gas_spec']):
            if spec[1] == Product:
                product_ind = spec[0]
        
        # Make sure the index has been found
        if product_ind == -1:
            print 'Product species not found'
        else:
            product_ind = product_ind + self.data.Species['n_surf']         # Adjust index to account for surface species                
        
        n_t_points = len(self.data.Specnum['t'])
        rate_traj = np.zeros(n_t_points)
        for t_point in range(n_t_points):
            if t_point == 0:
                rate_traj[t_point] = 0
            else:
                for i, elem_stoich in enumerate(self.data.Reactions['Nu']):
                    TOF_stoich = elem_stoich[product_ind]
                    r = self.data.Binary['propCounter'][t_point,i] / self.data.Specnum['t'][t_point]      # ergodic average
#                    r = (self.data.Binary['propCounter'][t_point,i] - self.data.Binary['propCounter'][t_point-1,i]) / (self.data.Specnum['t'][t_point] - self.data.Specnum['t'][t_point-1])      # non-ergodic average
                    rate_traj[t_point] = rate_traj[t_point] + TOF_stoich * r
        
        if rate_traj[-1] == 0:
            return False        
        
        t_vec = self.data.Specnum['t'][1::] / self.data.Specnum['t'][-1]
        rate_traj = (rate_traj[1::] - rate_traj[1]) / (rate_traj[-1] - rate_traj[1])
        
        n_back = int(n_t_points * (1-frac_sample))
        dydt = (rate_traj[-1] - rate_traj[n_back]) / (t_vec[-1] - t_vec[n_back]) 
        
        if show_graph:
            self.PlotOptions
            plt.figure()                 
            plt.plot(t_vec, rate_traj, color='r')
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.ylabel('observable',size=24)
            plt.xlabel('time',size=24)
            plt.show()        
        
        return np.abs(dydt) < d_cut