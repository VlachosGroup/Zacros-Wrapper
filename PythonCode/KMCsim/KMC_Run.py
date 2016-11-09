# -*- coding: utf-8 -*-
"""
Created on Thu Mar 03 14:54:26 2016

@author: robieta
"""

from IOdata import IOdata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mat

# For executable
import os
import sys
import subprocess
import copy
import matplotlib.animation as animation

class KMC_Run(IOdata):
    
    def __init__(self):
        
        super(KMC_Run, self).__init__()

        self.exe_file = ''
        self.anim = []          # animation object used for lattice movie
        self.props_avg = []
        self.rate_traj = []
        self.int_rate_traj = []
        
    def Run_sim(self):
        
        os.chdir(self.Path)
        
        print '--- Zacros run starting ---'
        subprocess.call([self.exe_file])
        print '--- Zacros run completed ---'

    def PlotOptions(self):
        mat.rcParams['mathtext.default'] = 'regular'
        mat.rcParams['text.latex.unicode'] = 'False'
        mat.rcParams['legend.numpoints'] = 1
        mat.rcParams['lines.linewidth'] = 2
        mat.rcParams['lines.markersize'] = 12

    def PlotSurfSpecVsTime(self, save = True):
        self.PlotOptions()
        plt.figure()
        
        for i in range (len(self.Species['surf_spec'])):
            plt.plot(self.Specnum['t'], self.Specnum['spec'][:,i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('spec. pop.',size=24)
#        plt.ylabel('coverage',size=24)
        plt.legend(self.Species['surf_spec'],loc=4,prop={'size':20},frameon=False)
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        
        if save:
            plt.savefig(self.Path + 'surf_spec_vs_time.png')
            plt.close()
        else:
            plt.show()
    
    def PlotGasSpecVsTime(self, save = True):
        self.PlotOptions()
        plt.figure()          
          
        for i in range (len(self.Species['gas_spec'])):
            ind = i + len(self.Species['surf_spec'])
            plt.plot(self.Specnum['t'], self.Specnum['spec'][:,ind])    
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('spec. pop.',size=24)
        plt.legend(self.Species['gas_spec'],loc=2,prop={'size':20},frameon=False)        
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)        
        
        if save:
            plt.savefig(self.Path + 'gas_spec_vs_time.png')
            plt.close()
        else:
            plt.show()
    
    def time_avg_props(self):
        
        delt = self.Specnum['t'][1::] - self.Specnum['t'][:-1:]
        props = self.Binary['propCounter'][1::,:] - self.Binary['propCounter'][:-1:,:]
        prop_shape = self.Binary['propCounter'].shape
        self.props_avg = np.zeros([prop_shape[0]-1, prop_shape[1]])
        
        for rxn_ind in range(prop_shape[1]):
            self.props_avg[:,rxn_ind] = props[:,rxn_ind] / delt
        
    def PlotPropsVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()
        
        self.time_avg_props()
        
        labels = []
        for i in range (self.props_avg.shape[1]):
            if np.max(np.abs(self.props_avg[:,i])) > 0:
#                plt.plot(self.Specnum['t'], self.Binary['prop'][:,i])
                plt.plot(self.Specnum['t'][1::], self.props_avg[:,i] )
                labels.append(self.Reactions['names'][i/2])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('props',size=24)
        plt.legend(labels,loc=2,prop={'size':20},frameon=False)
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        plt.show()
        
    def PlotIntPropsVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.Reactions['names'])):
            if np.max(np.abs(self.Binary['propCounter'][:,i])) > 0:
                plt.plot(self.Specnum['t'], self.Binary['propCounter'][:,i]) 
                labels.append(self.Reactions['names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('integral props',size=24)
        plt.legend(labels,loc=2,prop={'size':20},frameon=False)
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        plt.show()
    
    def PlotWVsTime(self):      # Helps analyze the sensitivty analysis
        self.PlotOptions
        plt.figure()            
            
        labels = []
        for i in range (len(self.Reactions['names'])):
            if np.max(np.abs(self.Binary['W_sen_anal'][:,i])) > 0:
                plt.plot(self.Specnum['t'], self.Binary['W_sen_anal'][:,i]) 
#                plt.plot(self.Specnum['t'], self.Procstat['events'][:,i] - self.Binary['propCounter'][:,i]) 
                labels.append(self.Reactions['names'][i])
        
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('W',size=24)
        plt.legend(labels, loc=2,prop={'size':20},frameon=False)
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        plt.show()       
    
    def PlotElemStepFreqs(self, save = True):
        self.PlotOptions
        plt.figure()        
        
        width = 0.2
        ind = 0
        yvals = []
        ylabels = []

        for i in range (self.Reactions['nrxns']):
            if self.Procstat['events'][-1,2*i] + self.Procstat['events'][-1,2*i+1] > 0:
                net_freq = abs(self.Procstat['events'][-1,2*i] - self.Procstat['events'][-1,2*i+1])
                if self.Procstat['events'][-1,2*i] > 0:              
                    plt.barh(ind-0.4, self.Procstat['events'][-1,2*i], width, color='r')
                if self.Procstat['events'][-1,2*i+1] > 0:
                    plt.barh(ind-0.6, self.Procstat['events'][-1,2*i+1], width, color='b')
                if net_freq > 0:
                    plt.barh(ind-0.8, net_freq, width, color='g')
                ylabels.append(self.Reactions['names'][i])
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
        
        if save:
            plt.savefig(self.Path + 'elem_step_freqs.png')
            plt.close()
        else:
            plt.show()
    
    def LatticeMovie(self):       # Need to complete this function by plotting adsorbates from the history file data

        self.KMC_lat.Read_lattice_output(self.Path + 'lattice_output.txt')
        cart_coords = self.KMC_lat.cart_coords
        lat = self.KMC_lat.lattice_matrix
#        self.KMC_lat.PlotLattice()
        border = np.dot(np.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]), lat)
        
        self.PlotOptions

        fig = plt.figure()
        ax = plt.axes(xlim=[np.min(border[:,0]), np.max(border[:,0])], ylim=[np.min(border[:,1]), np.max(border[:,1])])
        plt.axis('equal')        
        
        # Plot cell border 
        plt.plot(border[:,0], border[:,1], '--k', linewidth = 4)
        
        # Plot neighbors
        cutoff = 3.0
        for pair in self.KMC_lat.neighbor_list:
                p1 = np.array([cart_coords[pair[0]-1,0], cart_coords[pair[0]-1,1]])
                p2 = np.array([cart_coords[pair[1]-1,0], cart_coords[pair[1]-1,1]])
                if np.linalg.norm(p2 - p1) < cutoff:
                    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], '-k', linewidth = 1)
        
        # Plot sites            
        plt.plot(cart_coords[:,0], cart_coords[:,1], 'bo', markersize = 15)
            
        color_list = ['g','r','c','m','y','k']
        line_list = []
        for ind in range(self.Species['n_surf']):    
            line, = ax.plot([], [], 's' + color_list[ind % len(color_list)], markersize = 10, label=self.Species['surf_spec'][ind])
            line_list.append(line)
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('x-coord (ang)',size=24)
        plt.ylabel('y-coord (ang)',size=24)        
        plt.legend(handles = line_list, frameon=False)
        
        # initialization function: plot the background of each frame
        def init():
            for line in line_list:
                line.set_data([], [])
            return line_list
        
        # animation function.  This is called sequentially
        def animate(i):
            
            snap = self.History[i]            
            
            for ind in range(self.Species['n_surf']):
                
                # Find all coordinates with species ind occupying it
                x_list = []
                y_list = []
                for site_ind in range(self.n_sites):
                    if snap[site_ind,2] == ind + 1:
                        x_list.append(cart_coords[site_ind,0])
                        y_list.append(cart_coords[site_ind,1])
        
                x = np.array(x_list)
                y = np.array(y_list)                
                
                line = line_list[ind]
                line.set_data(x,y)
                
            return line_list
        
        # call the animator.  blit=True means only re-draw the parts that have changed.
        self.anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=self.n_snapshots, interval=40, blit=True, repeat_delay = 500)
        
        plt.show()
    
    def ComputeTOF(self,Product):                       # return TOF and TOF error
        
        # Find the index of the product species
        product_ind = -1       
        for i in enumerate(self.Species['gas_spec']):
            if i[1] == Product:
                product_ind = i[0]
        
        # Make sure the index has been found
        if product_ind == -1:
            print 'Product species not found'
        else:
            product_ind = product_ind + self.Species['n_surf']         # Adjust index to account for surface species
        
        
        nRxns = len(self.Reactions['Nu'])        
        TOF_contributions = [0 for i in range(nRxns)]              # number of product molecules produced in each reaction        
        for i, elem_stoich in enumerate(self.Reactions['Nu']):
            TOF_stoich = elem_stoich[product_ind]
            r = self.Binary['propCounter'][-1,i] / self.Specnum['t'][-1]      # ergodic average
#            r = self.Binary['prop'][-1,i]                                           # non-ergodic average
            TOF_contributions[i] = TOF_stoich * r         
               
        TOF = np.sum(TOF_contributions)
        TOF_fracs = TOF_contributions / TOF             # will need this for sensitivity analysis
#        return TOF
        return {'TOF': TOF, 'TOF_fracs': TOF_fracs}    
    
    def AdjustPreExponentials(self, delta_sdf):
        
        rxn_ind = 0
        for rxn_type in self.Reactions['Input']:
            for variant in rxn_type['variant']:
                variant['pre_expon'] = variant['pre_expon'] * delta_sdf[rxn_ind]
                self.scaledown_factors[rxn_ind] = self.scaledown_factors[rxn_ind] * delta_sdf[rxn_ind]
                rxn_ind += 1    
    
    def CalcRateTraj(self, Product):
        
        try:
            product_ind = self.Species['n_surf'] + self.Species['gas_spec'].index(Product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Gas species not found.')
        
        n_t_points = len(self.Specnum['t'])
        self.rate_traj = np.zeros(n_t_points)
        self.int_rate_traj = np.zeros(n_t_points)
        
        for t_point in range(n_t_points):
            for i, elem_stoich in enumerate(self.Reactions['Nu']):
                
                TOF_stoich = elem_stoich[product_ind]
                r_int = self.Binary['propCounter'][t_point,i]
                self.int_rate_traj[t_point] = self.int_rate_traj[t_point] + TOF_stoich * r_int
                
                if t_point == 0:
                    self.rate_traj[t_point] = 0
                else:                    
                    r = (self.Binary['propCounter'][t_point,i] - self.Binary['propCounter'][t_point-1,i]) / (self.Specnum['t'][t_point] - self.Specnum['t'][t_point-1])      # ergodic average
                    self.rate_traj[t_point] = self.rate_traj[t_point] + TOF_stoich * r
                        
        self.rate_traj = self.rate_traj[1::]
    
    def PlotRateVsTime(self):
    
        self.PlotOptions()
        
        plt.figure()
        plt.plot(self.Specnum['t'][1::], self.rate_traj)
        plt.xticks(size=20)
        plt.yticks(size=20)
        plt.xlabel('time (s)',size=24)
        plt.ylabel('rate',size=24)
        plt.ylim([0, plt.ylim()[1]])
        ax = plt.subplot(111)
        pos = [0.2, 0.15, 0.7, 0.8]
        ax.set_position(pos)
        plt.show()    
    
    def time_search(self, t):
        
        if t > self.Specnum['t'][-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind = 0
        while self.Specnum['t'][ind] < t:
            ind += 1
            
        return ind
    
    def avg_in_window(self, data, limits):
        high_ind = self.time_search(limits[1])
        low_ind = self.time_search(limits[0])
        return (data[high_ind] - data[low_ind]) / (self.Specnum['t'][high_ind] - self.Specnum['t'][low_ind])
    
    def CheckSteadyState(self, Product, frac_sample = 0.2, cut = 0.05, show_graph = False):
           
        self.CalcRateTraj(Product)
        
        if self.rate_traj[-1] == 0:          # Simulation is not in steady-state if rate of production of product species is 0
            return False        

        batch_1 = self.avg_in_window(self.int_rate_traj, [0.50 * self.Specnum['t'][-1], 0.75 * self.Specnum['t'][-1]])
        batch_2 = self.avg_in_window(self.int_rate_traj, [0.75 * self.Specnum['t'][-1], self.Specnum['t'][-1]])
        
        frac_change = (batch_2 - batch_1) / batch_2
        
        print '{0:.3f}'.format(frac_change) + ' change in last 50% of window'
        
        if show_graph:
            self.PlotOptions
            plt.figure()                 
            plt.plot(self.Specnum['t'][1::], self.rate_traj, color='r')
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.ylabel('rate',size=24)
            plt.xlabel('time (s)',size=24)
            plt.ylim([0, plt.ylim()[1]])
            plt.show()
            
        return np.abs(frac_change) < cut
        
    @staticmethod
    def time_sandwich(run1, run2):
        
        sandwich = copy.deepcopy(run1)
        sandwich.Performance['t_final'] = run1.Performance['t_final'] + run2.Performance['t_final']
        sandwich.Performance['events_occurred'] = run1.Performance['events_occurred'] + run2.Performance['events_occurred']
        sandwich.Performance['CPU_time'] = run1.Performance['CPU_time'] + run2.Performance['CPU_time']
        
        sandwich.Specnum['t'] = np.concatenate([run1.Specnum['t'], run2.Specnum['t'][1::] + run1.Performance['t_final'] * np.ones( len(run2.Specnum['t'])-1 )])
        sandwich.Specnum['spec'] = np.vstack([run1.Specnum['spec'], run2.Specnum['spec'][1::,:] ])
#        print run1.Procstat['events']
#        print run2.Procstat['events'][1::,:]
#        print run1.Procstat['events'][-1,:]
#        print np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), run1.Procstat['events'][-1,:])
        sandwich.Procstat['events'] = np.vstack( [run1.Procstat['events'], run2.Procstat['events'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Procstat['events'][-1,:]] ) ] )
        sandwich.Binary['propCounter'] = np.vstack( [run1.Binary['propCounter'], run2.Binary['propCounter'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Binary['propCounter'][-1,:]]  ) ] )
        
        return sandwich