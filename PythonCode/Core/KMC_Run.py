from IOdata import IOdata
import numpy as np
import matplotlib as mat
mat.use('Agg')
import matplotlib.pyplot as plt

# For executable
import os
import sys
import subprocess
import copy
import matplotlib.animation as animation
from Helper import FileIO

from scipy.optimize import curve_fit

# Use this functional form to regress the gas profiles
def func(t, rss, tau):

    t = np.asarray(t)
    y = np.zeros(t.shape)
    y += (t / tau >= 10) * rss * ( t - tau )
    y += (t / tau < 10) * rss * ( t - tau * ( 1 - np.exp( - t / tau) ) )

    return y

class KMC_Run(IOdata):
    
    def __init__(self):
        
        super(KMC_Run, self).__init__()

        self.exe_file = ''
        self.anim = []          # animation object used for lattice movie
        self.props_avg = []
        self.rate_traj = []
        self.int_rate_traj = []
        self.gas_stoich = []        # stoichiometry of gas-phase reaction
        self.net_rxn = []        # net mass, conserved at steady-state
        
    def Run_sim(self):
        
        os.chdir(self.Path)
        
        if self.exe_file == '':
            raise Exception('Zacros exe file not specified.')        
        
        try:
            print '--- Zacros run starting ---'
            subprocess.call([self.exe_file])
            print '--- Zacros run completed ---'
        except:
            raise Exception('Zacros run failed.')
    
    def ComputeTOF(self, Product, win = [0.0, 1.0]):                       # return TOF and TOF error
        
        # Find the index of the product species
        try:
            product_ind = self.Species['n_surf'] + self.Species['gas_spec'].index(Product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species not found.')
        
        # Find indices for the beginning and end of the time window
        start_t = win[0] * self.Specnum['t'][-1]
        end_t = win[1] * self.Specnum['t'][-1]
        start_ind = self.time_search(start_t)
        end_ind = self.time_search(end_t)
        del_t = self.Specnum['t'][end_t] - self.Specnum['t'][start_ind]          
        
        # Compute rates based on propensities and stoichiometries
        nRxns = len(self.Reactions['Nu'])
        TOF_contributions_inst = np.zeros(nRxns)
        TOF_contributions_erg = np.zeros(nRxns)
        for i, elem_stoich in enumerate(self.Reactions['Nu']):
            TOF_stoich = elem_stoich[product_ind]
            r_inst = ( self.Binary['propCounter'][end_ind,i] - self.Binary['propCounter'][end_ind-1,i] ) / ( self.Specnum['t'][end_ind] - self.Specnum['t'][end_ind-1] )
            r_erg = ( self.Binary['propCounter'][end_ind,i] - self.Binary['propCounter'][start_ind,i] ) / del_t      # ergodic average
            TOF_contributions_inst[i] = TOF_stoich * r_inst            
            TOF_contributions_erg[i] = TOF_stoich * r_erg         
        
        # Total the rates
        TOF_inst = np.sum(TOF_contributions_inst)
        TOF_fracs_inst = TOF_contributions_inst / TOF_inst
        TOF_erg = np.sum(TOF_contributions_erg)
        TOF_fracs_erg = TOF_contributions_erg / TOF_erg

        return {'TOF_inst': TOF_inst, 'TOF_fracs_inst': TOF_fracs_inst, 'TOF_erg': TOF_erg, 'TOF_fracs_erg': TOF_fracs_erg}    
    
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
            raise Exception('Product species not found.')
        
        n_t_points = len(self.Specnum['t'])
        self.rate_traj = np.zeros(n_t_points)
        self.int_rate_traj = np.zeros(n_t_points)
        
        for t_point in range(n_t_points):
            for i, elem_stoich in enumerate(self.Reactions['Nu']):
                
                TOF_stoich = elem_stoich[product_ind]
                
                if t_point == 0:
                    self.rate_traj[t_point] = 0
                    self.int_rate_traj[t_point] = 0
                else:                    
                    r = (self.Binary['propCounter'][t_point,i] - self.Binary['propCounter'][t_point-1,i]) / (self.Specnum['t'][t_point] - self.Specnum['t'][t_point-1])      # averaged in an interval
                    r_int = self.Binary['propCounter'][t_point,i] / self.Specnum['t'][t_point]      # ergodic average              
                    self.rate_traj[t_point] = self.rate_traj[t_point] + TOF_stoich * r
                    self.int_rate_traj[t_point] = self.int_rate_traj[t_point] + TOF_stoich * r_int
                        
        self.rate_traj = self.rate_traj[1::]
        self.int_rate_traj = self.int_rate_traj[1::]
    
    def time_search(self, t):           # Convert KMC time to index of time series
        
        if t > self.Specnum['t'][-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind = 0
        while self.Specnum['t'][ind] < t:
            ind += 1
            
        return ind
        
    def fraction_search(self, frac):       # Convert percentage of time interval to index of time series
        
        t = frac * self.Specnum['t'][-1]
        return self.time_search(t)
        
    def avg_in_window(self, data, limits):
        high_ind = self.time_search(limits[1])
        low_ind = self.time_search(limits[0])
        return (data[high_ind] - data[low_ind]) / (self.Specnum['t'][high_ind] - self.Specnum['t'][low_ind])
    
    def CheckSteadyState(self, Product, frac_sample = 0.2, cut = 0.05, show_graph = False):
           
        self.CalcRateTraj(Product)
        
        if self.rate_traj[-1] == 0:          # Simulation is not in steady-state if rate of production of product species is 0
            return False        

        batch_1 = self.avg_in_window(self.int_rate_traj, [0, 0.5 * self.Specnum['t'][-1]])
        batch_2 = self.avg_in_window(self.int_rate_traj, [0.5 * self.Specnum['t'][-1], self.Specnum['t'][-1]])
        
        frac_change = (batch_2 - batch_1) / batch_2
        
        print '{0:.3f}'.format(frac_change) + ' change in last 50% of window'
            
        return np.abs(frac_change) < cut
        
    @staticmethod
    def time_sandwich(run1, run2):
        
        sandwich = copy.deepcopy(run1)
        sandwich.Performance['t_final'] = run1.Performance['t_final'] + run2.Performance['t_final']
        sandwich.Performance['events_occurred'] = run1.Performance['events_occurred'] + run2.Performance['events_occurred']
        sandwich.Performance['CPU_time'] = run1.Performance['CPU_time'] + run2.Performance['CPU_time']
        
        sandwich.Specnum['t'] = np.concatenate([run1.Specnum['t'], run2.Specnum['t'][1::] + run1.Specnum['t'][-1] * np.ones( len(run2.Specnum['t'])-1 )])
        
        n_surf_specs = len(run2.Species['surf_spec'])
        run2.Specnum['spec'][1::, n_surf_specs : ] = run2.Specnum['spec'][1::, n_surf_specs : ] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Specnum['spec'][-1, n_surf_specs : ]] )        
        sandwich.Specnum['spec'] = np.vstack([run1.Specnum['spec'], run2.Specnum['spec'][1::,:] ])
        sandwich.Procstat['events'] = np.vstack( [run1.Procstat['events'], run2.Procstat['events'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Procstat['events'][-1,:]] ) ] )
        sandwich.Binary['propCounter'] = np.vstack( [run1.Binary['propCounter'], run2.Binary['propCounter'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Binary['propCounter'][-1,:]]  ) ] )
        sandwich.Binary['W_sen_anal']  = np.vstack( [run1.Binary['W_sen_anal'], run2.Binary['W_sen_anal'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Binary['W_sen_anal'][-1,:]]  ) ] )      
        
        sandwich.History = [ run1.History[0], run2.History[-1] ]     
        
        return sandwich
    
    def time_avg_props(self):
        
        delt = self.Specnum['t'][1::] - self.Specnum['t'][:-1:]
        props = self.Binary['propCounter'][1::,:] - self.Binary['propCounter'][:-1:,:]
        prop_shape = self.Binary['propCounter'].shape
        self.props_avg = np.zeros([prop_shape[0]-1, prop_shape[1]])
        
        for rxn_ind in range(prop_shape[1]):
            self.props_avg[:,rxn_ind] = props[:,rxn_ind] / delt
            
    def calc_net_rxn(self):
        
        gas_mat = self.Specnum['spec'][:, len(self.Species['surf_spec']) : ]        # gas species only
        self.net_rxn = np.dot(gas_mat, np.abs(self.gas_stoich))
    
    def CheckNetRxnConvergence(self, m = 0.04):
        
        half_ind = self.fraction_search(0.5)
        diff1 = np.abs( (self.net_rxn[half_ind] - self.net_rxn[-1]) / (self.net_rxn[-1] - self.net_rxn[0]) )

        return diff1 < m
        
    def CheckProductConvergence(self, Product, dimless_time_limit = 12.0):
        
        try:
            product_ind = self.Species['n_surf'] + self.Species['gas_spec'].index(Product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species not found.')
        
        tau_fit = self.FitProductProfile(Product)
		
        if self.Specnum['spec'][:,product_ind][-1] <= 100:
            return False
        else:
            return self.Specnum['t'][-1] / tau_fit > dimless_time_limit

    def FitProductProfile(self, Product, graph_compare = True):

        try:
            product_ind = self.Species['n_surf'] + self.Species['gas_spec'].index(Product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species not found.')
    
        rss_guess = self.Specnum['spec'][:,product_ind][-1] / self.Specnum['t'][-1]
        tau_guess = self.Specnum['t'][-1] / 5.0
    
        popt, pcov = curve_fit(func, self.Specnum['t'], self.Specnum['spec'][:,product_ind], p0 = [rss_guess, tau_guess], bounds=[ np.array([-np.inf, 0]), np.array([np.inf, np.inf]) ])

        if graph_compare:

            spec_fit = func(self.Specnum['t'], popt[0], popt[1])

            mat.rcParams['mathtext.default'] = 'regular'
            mat.rcParams['text.latex.unicode'] = 'False'
            mat.rcParams['legend.numpoints'] = 1
            mat.rcParams['lines.linewidth'] = 2
            mat.rcParams['lines.markersize'] = 12
    
            plt.figure()
    
            plt.plot(self.Specnum['t'], self.Specnum['spec'][:,product_ind])
            plt.plot(self.Specnum['t'], spec_fit)
    
    
            plt.xticks(size=20)
            plt.yticks(size=20)
            plt.xlabel('Time (s)', size=24)
            plt.ylabel(Product + ' count', size=24)
    
            plt.legend(['data', 'fit'], loc=4, prop={'size':20}, frameon=False)
    
            ax = plt.subplot(111)
            ax.set_position([0.2, 0.15, 0.7, 0.8])
    
            plt.savefig('fit_compare.png')
            
        return popt[1]

    ''' ==================================== Plotting methods ==================================== '''
    
    def PlotSurfSpecVsTime(self):
        
        time_vecs = []
        surf_spec_vecs = []
        for i in range (len(self.Species['surf_spec'])):
            time_vecs.append(self.Specnum['t'])
            surf_spec_vecs.append(self.Specnum['spec'][:,i])
        
        FileIO.PlotTrajectory(time_vecs, surf_spec_vecs, xlab = 'Time (s)', ylab = 'spec. pop.', series_labels = self.Species['surf_spec'], fname = os.path.join(self.Path, 'surf_spec_vs_time.png'))
    
    def PlotGasSpecVsTime(self):
        
        time_vecs = []
        gas_spec_vecs = []
        for i in range (len(self.Species['gas_spec'])):
            time_vecs.append(self.Specnum['t'])
            gas_spec_vecs.append(self.Specnum['spec'][:,i + len(self.Species['surf_spec']) ])
        
        FileIO.PlotTrajectory(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'spec. pop.', series_labels = self.Species['gas_spec'], fname = os.path.join(self.Path, 'gas_spec_vs_time.png'))
        
    def PlotNetGasRxnVsTime(self):
        FileIO.PlotTrajectory([self.Specnum['t']], [self.net_rxn], xlab = 'Time (s)', ylab = 'spec. pop.', fname = os.path.join(self.Path, 'net_rxn_vs_time.png'))
        
    def PlotPropsVsTime(self):
        
        self.time_avg_props()
        time_vecs = []
        gas_spec_vecs = []
        labels = []
        for i in range (self.props_avg.shape[1]):
            time_vecs.append(self.Specnum['t'][1::])
            gas_spec_vecs.append(self.props_avg[:,i])
            labels.append(self.Reactions['names'][i/2])
        
        FileIO.PlotTrajectory(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'props (1/s)', series_labels = labels, fname = os.path.join(self.Path, 'props_vs_time.png'))
        
    def PlotIntPropsVsTime(self, save = True):      # Helps analyze the sensitivty analysis
        
        self.time_avg_props()
        time_vecs = []
        gas_spec_vecs = []
        labels = []
        for i in range (self.props_avg.shape[1]):
            time_vecs.append(self.Specnum['t'])
            gas_spec_vecs.append(self.Binary['propCounter'][:,i])
            labels.append(self.Reactions['names'][i/2])
        
        FileIO.PlotTrajectory(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'int props', series_labels = labels, fname = os.path.join(self.Path + 'int_props_vs_time.png'))
        
        FileIO.PlotOptions
        plt.figure()
    
    def PlotWVsTime(self):      # Helps analyze the sensitivty analysis
    
        time_vecs = []
        W_vecs = []
        labels = []
        for i in range (len(self.Reactions['names'])):
            if np.max(np.abs(self.Binary['W_sen_anal'][:,i])) > 0:
                time_vecs.append(self.Specnum['t'])
                W_vecs.append( self.Binary['W_sen_anal'][:,i])
                labels.append(self.Reactions['names'][i])
        
        FileIO.PlotTrajectory(time_vecs, W_vecs, xlab = 'Time (s)', ylab = 'traj. deriv.', series_labels = labels, fname = os.path.join(self.Path + 'traj_deriv_vs_time.png'))
    
    def PlotElemStepFreqs(self, window = [0.0, 1.0], time_norm = False, site_norm = 1.0):
        
        start_ind = self.fraction_search(window[0])
        end_ind = self.fraction_search(window[1])
        event_freqs = ( self.Procstat['events'][end_ind,:] - self.Procstat['events'][start_ind,:] ) / site_norm
        if time_norm:
            event_freqs = event_freqs / ( self.Specnum['t'][end_ind] - self.Specnum['t'][start_ind] )
        
        FileIO.PlotOptions
        plt.figure()        
        
        width = 0.2
        ind = 0
        yvals = []
        ylabels = []
        bar_vals = []
        for i in range (self.Reactions['nrxns']):
            if event_freqs[2*i] + event_freqs[2*i+1] > 0:
                net_freq = abs(event_freqs[2*i] - event_freqs[2*i+1])
                if event_freqs[2*i] > 0:              
                    plt.barh(ind-0.4, event_freqs[2*i], width, color='r', log = True)
                    bar_vals.append(event_freqs[2*i])
                if event_freqs[2*i+1] > 0:
                    plt.barh(ind-0.6, event_freqs[2*i+1], width, color='b', log = True)
                    bar_vals.append(event_freqs[2*i+1])
                if net_freq > 0:
                    plt.barh(ind-0.8, net_freq, width, color='g', log = True)
                    bar_vals.append(net_freq)
                ylabels.append(self.Reactions['names'][i])
                yvals.append(ind-0.6)
                ind = ind - 1

        bar_vals = np.array(bar_vals)
        log_bar_vals = np.log10(bar_vals)
        xmin = 10**np.floor(np.min(log_bar_vals))
        xmax = 10**np.ceil(np.max(log_bar_vals))
                
        plt.xticks(size=20)
        plt.yticks(size=10)
        plt.xlabel('Frequency',size=24)
        plt.yticks(yvals, ylabels)
        plt.legend(['fwd','bwd','net'],loc=4,prop={'size':20},frameon=False)
        plt.xlim([xmin, xmax])
        ax = plt.subplot(111)        
        pos = [0.30, 0.15, 0.65, 0.8]
        ax.set_position(pos)
        
        plt.savefig(os.path.join(self.Path, 'elem_step_freqs.png'))
        plt.close()


    def PrintElemStepFreqs(self, window = [0.0, 1.0], time_norm = False, site_norm = 1.0):      # into an output file
        
        start_ind = self.fraction_search(window[0])
        end_ind = self.fraction_search(window[1])
        event_freqs = ( self.Procstat['events'][end_ind,:] - self.Procstat['events'][start_ind,:] ) / site_norm
        if time_norm:
            event_freqs = event_freqs / ( self.Specnum['t'][end_ind] - self.Specnum['t'][start_ind] )
                
            
        with open(os.path.join(self.Path, 'rxn_freqs.txt'), 'w') as txt:   
            txt.write('----- Elementary reaction frequencies -----\n')
            txt.write('Reaction Name \t Forward \t Reverse \t Net \n')
            for i in range (self.Reactions['nrxns']):
                net_freq = abs(event_freqs[2*i] - event_freqs[2*i+1])
                txt.write(self.Reactions['names'][i] + '\t')
                txt.write('{0:.3E} \t'.format(event_freqs[2*i]))
                txt.write('{0:.3E} \t'.format(event_freqs[2*i+1]))
                txt.write('{0:.3E} \n'.format(net_freq))


    def PlotRateVsTime(self):
        FileIO.PlotTrajectory([self.Specnum['t'][1::]], [self.rate_traj], xlab = 'Time (s)', ylab = 'rate (1/s)', fname = os.path.join(self.Path + 'rate_vs_time.png'))
        FileIO.PlotTrajectory([self.Specnum['t'][1::]], [self.int_rate_traj], xlab = 'Time (s)', ylab = 'rate (1/s)', fname = os.path.join(self.Path + 'rate_erg_vs_time.png'))

    def LatticeMovie(self):       # Need to complete this function by plotting adsorbates from the history file data

        self.KMC_lat.Read_lattice_output(os.path.join(self.Path, 'lattice_output.txt'))
        cart_coords = self.KMC_lat.cart_coords
        lat = self.KMC_lat.lattice_matrix
#        self.KMC_lat.PlotLattice()
        border = np.dot(np.array([[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]), lat)
        
        FileIO.PlotOptions

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