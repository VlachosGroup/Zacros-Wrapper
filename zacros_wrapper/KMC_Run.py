from IOdata import IOdata
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# For executable
import os
import sys
import subprocess
import copy

from Helper import *
import scipy

class kmc_traj(IOdata):
    
    '''
    Handles a single Zacros trajectory. Inherits from IOdata.
    '''
    
    def __init__(self):
        
        '''
        Initializes additional class variables
        '''
        
        super(kmc_traj, self).__init__()

        self.exe_file = None
        self.props_avg = None
        self.rate_traj = None
        self.int_rate_traj = None
        self.gas_stoich = None        # stoichiometry of gas-phase reaction
        
    def Run_sim(self):
        
        '''
        Call the Zacros executable and run the simulation
        '''
        
        os.chdir(self.Path)
        
        if self.exe_file is None:
            raise Exception('Zacros executable not specified.')        
        
        try:
            print '--- Zacros run starting ---'
            subprocess.call([self.exe_file])
            print '--- Zacros run completed ---'
        except:
            raise Exception('Zacros run failed.')
    
    def ComputeTOF(self, Product, win = [0.0, 1.0]):                       # return TOF and TOF error
        
        '''
        Compute the turnover frequency over a given window
        '''
        
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
        
        '''
        Adjust the pre-exponential ratios of all elementary reactions
        '''
        
        rxn_ind = 0
        for rxn_type in self.Reactions['Input']:
            for variant in rxn_type['variant']:
                variant['pre_expon'] = variant['pre_expon'] * delta_sdf[rxn_ind]
                self.scaledown_factors[rxn_ind] = self.scaledown_factors[rxn_ind] * delta_sdf[rxn_ind]
                rxn_ind += 1    
    
    def CalcRateTraj(self, Product):
        
        '''
        Compute instantaneous rates
        '''
        
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
    
    def time_search(self, t):
        
        '''
        Given a time, look up the index of the time vector nearest to that time
        '''
        
        if t > self.Specnum['t'][-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind = 0
        while self.Specnum['t'][ind] < t:
            ind += 1
            
        return ind
        
        
    def avg_in_window(self, data, limits):
    
        '''
        Ergodically average data within a time window
        '''
    
        high_ind = self.time_search(limits[1])
        low_ind = self.time_search(limits[0])
        return (data[high_ind] - data[low_ind]) / (self.Specnum['t'][high_ind] - self.Specnum['t'][low_ind])
    
        
    @staticmethod
    def time_sandwich(run1, run2):
        
        '''
        Take two trajectories and append them
        '''
        
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
        
        '''
        ???
        '''
        
        delt = self.Specnum['t'][1::] - self.Specnum['t'][:-1:]
        props = self.Binary['propCounter'][1::,:] - self.Binary['propCounter'][:-1:,:]
        prop_shape = self.Binary['propCounter'].shape
        self.props_avg = np.zeros([prop_shape[0]-1, prop_shape[1]])
        
        for rxn_ind in range(prop_shape[1]):
            self.props_avg[:,rxn_ind] = props[:,rxn_ind] / delt 
    
        
    def CheckGasConvergence(self, gas_stoich, r_sqr_cut = 0.98, rate_perc = 0.03):
    
        '''
        Check whether the simulation is at steady state using gas species profiles
        '''
    
        window = [0.0, 1.0]
        start_ind = self.fraction_search(window[0])
        end_ind = self.fraction_search(window[1])
        rates = []
        r_sqr_vals = []
    
        for gas_spec_ind in range (len(self.Species['gas_spec'])):
        
            x = self.Specnum['t'][start_ind:end_ind]
            y = self.Specnum['spec'][start_ind:end_ind, gas_spec_ind + len(self.Species['surf_spec']) ]

            (a_s, b_s, r, tt, stderr) = scipy.stats.linregress(x, y)
            print('Linear regression using stats.linregress')
            print('\nregression: a=%.2f b=%.2f, std error= %.3f' % (a_s, b_s, stderr))
            print('r^2 value: %.2f' % (r**2))
            
            mat.rcParams['mathtext.default'] = 'regular'
            mat.rcParams['text.latex.unicode'] = 'False'
            mat.rcParams['legend.numpoints'] = 1
            mat.rcParams['lines.linewidth'] = 2
            mat.rcParams['lines.markersize'] = 12
                    
            plt.figure()
            
            plt.plot(self.Specnum['t'], self.Specnum['spec'][:, gas_spec_ind + len(self.Species['surf_spec']) ])
            plt.plot(self.Specnum['t'], self.Specnum['t'] * a_s + b_s)
            #plt.plot(self.Specnum['t'], self.Specnum['spec'][:, gas_spec_ind + len(self.Species['surf_spec']) ] - (self.Specnum['t'] * a_s + b_s))
            
            plt.xlabel('Time (s)', size=24)
            plt.ylabel('Gas population', size=24)
            plt.legend(['data', 'fit'], loc=1, prop={'size':20}, frameon=False)
            
            ax = plt.subplot(111)
            ax.set_position([0.2, 0.15, 0.7, 0.8])
            
            plt.savefig(os.path.join(self.Path, self.Species['gas_spec'][gas_spec_ind] + '_fit.png'))
            plt.close()
            
            if gas_stoich[gas_spec_ind] != 0:
                rates.append(a_s / gas_stoich[gas_spec_ind])
                r_sqr_vals.append(r**2)
                
        print rates
        print r_sqr_vals    
            

    ''' ==================================== Plotting methods ==================================== '''
    
    def PlotSurfSpecVsTime(self, site_norm = 1):
        
        '''
        Plot surface species profiles versus time
        '''
        
        if site_norm == 1:
            ylabel = 'Species count'
        else:
            ylabel = 'Coverage'
        
        time_vecs = []
        surf_spec_vecs = []
        for i in range (len(self.Species['surf_spec'])):
            time_vecs.append(self.Specnum['t'])
            surf_spec_vecs.append(self.Specnum['spec'][:,i] / float(site_norm))
        
        PlotTimeSeries(time_vecs, surf_spec_vecs, xlab = 'Time (s)', ylab = ylabel, series_labels = self.Species['surf_spec'], fname = os.path.join(self.Path, 'surf_spec_vs_time.png'))
    
    def PlotGasSpecVsTime(self):
        
        '''
        Plot gas phase species profiles versus time
        '''
        
        time_vecs = []
        gas_spec_vecs = []
        for i in range (len(self.Species['gas_spec'])):
            time_vecs.append(self.Specnum['t'])
            gas_spec_vecs.append(self.Specnum['spec'][:,i + len(self.Species['surf_spec']) ])
        
        PlotTimeSeries(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'spec. pop.', series_labels = self.Species['gas_spec'], fname = os.path.join(self.Path, 'gas_spec_vs_time.png'))
        
        
    def PlotPropsVsTime(self):
        
        '''
        Plot elementary step propensities versus time
        '''
        
        self.time_avg_props()
        time_vecs = []
        gas_spec_vecs = []
        labels = []
        for i in range (self.props_avg.shape[1]):
            time_vecs.append(self.Specnum['t'][1::])
            gas_spec_vecs.append(self.props_avg[:,i])
            labels.append(self.Reactions['names'][i/2])
        
        PlotTimeSeries(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'props (1/s)', series_labels = labels, fname = os.path.join(self.Path, 'props_vs_time.png'))
        
    def PlotIntPropsVsTime(self, save = True):      # Helps analyze the sensitivty analysis
        
        '''
        Plot integral propensities versus time
        '''
        
        self.time_avg_props()
        time_vecs = []
        gas_spec_vecs = []
        labels = []
        for i in range (self.props_avg.shape[1]):
            time_vecs.append(self.Specnum['t'])
            gas_spec_vecs.append(self.Binary['propCounter'][:,i])
            labels.append(self.Reactions['names'][i/2])
        
        PlotTimeSeries(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'int props', series_labels = labels, fname = os.path.join(self.Path + 'int_props_vs_time.png'))
        
        PlotOptions
        plt.figure()
    
    def PlotWVsTime(self):      # Helps analyze the sensitivty analysis
    
        '''
        Plot trajectory derivatives versus time
        '''
    
        time_vecs = []
        W_vecs = []
        labels = []
        for i in range (len(self.Reactions['names'])):
            if np.max(np.abs(self.Binary['W_sen_anal'][:,i])) > 0:
                time_vecs.append(self.Specnum['t'])
                W_vecs.append( self.Binary['W_sen_anal'][:,i])
                labels.append(self.Reactions['names'][i])
        
        PlotTimeSeries(time_vecs, W_vecs, xlab = 'Time (s)', ylab = 'traj. deriv.', series_labels = labels, fname = os.path.join(self.Path + 'traj_deriv_vs_time.png'))
    
    def PlotElemStepFreqs(self, window = [0.0, 1.0], time_norm = False, site_norm = 1):
        
        '''
        Plot a bar graph of elementary step frequencies versus time
        '''
        
        start_ind = self.fraction_search(window[0])
        end_ind = self.fraction_search(window[1])
        event_freqs = ( self.Procstat['events'][end_ind,:] - self.Procstat['events'][start_ind,:] ) / float(site_norm)
        if time_norm:
            event_freqs = event_freqs / ( self.Specnum['t'][end_ind] - self.Specnum['t'][start_ind] )
        
        PlotOptions
        plt.figure()
        
        width = 0.2
        ind = 0
        yvals = []
        ylabels = []
        bar_vals = []
        
        store_ind = 0       # index of where the data is stored
        for i in range (self.Reactions['nrxns']):
        
            fwd_rate = event_freqs[store_ind]
            store_ind += 1
            bwd_rate = 0
            
            if self.Reactions['is_reversible'][i]:
                bwd_rate = event_freqs[store_ind]
                store_ind += 1
        
            if fwd_rate + bwd_rate > 0:
            
                net_freq = abs(fwd_rate - bwd_rate)
                
                if fwd_rate > 0:              
                    plt.barh(ind-0.4, fwd_rate, width, color='r', log = True)
                    bar_vals.append(fwd_rate)
                if bwd_rate > 0:
                    plt.barh(ind-0.6, bwd_rate, width, color='b', log = True)
                    bar_vals.append(bwd_rate)
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
        
        '''
        Write elementary step frequencies into an output file
        '''
        
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
    
        '''
        Plot instantaneous rate versus time
        '''
    
        PlotTimeSeries([self.Specnum['t'][1::]], [self.rate_traj], xlab = 'Time (s)', ylab = 'rate (1/s)', fname = os.path.join(self.Path + 'rate_vs_time.png'))
        PlotTimeSeries([self.Specnum['t'][1::]], [self.int_rate_traj], xlab = 'Time (s)', ylab = 'rate (1/s)', fname = os.path.join(self.Path + 'rate_erg_vs_time.png'))

        
    def PlotLattice(self):
    
        '''
        Plot the lattice
        '''
    
        plt = self.KMC_lat.PlotLattice()
        plt.savefig(os.path.join(self.Path, 'lattice.png'))
        plt.close()
        
    def LatticeMovie(self, include_neighbor_lines = False, spec_color_list = ['b', 'g','r','c','m','y','k']):       # Need to complete this function by plotting adsorbates from the history file data

        '''
        Create a .png file with a picture of a lattice for every snapshot in history_output.txt
        '''
    
        cart_coords = self.KMC_lat.cart_coords
        spec_label_list = self.Species['surf_spec']
        
        frame_fldr = os.path.join(self.Path, 'lattice_frames')
        if not os.path.exists( frame_fldr ):
                os.makedirs( frame_fldr )
                
        print str(self.n_snapshots) + ' total snapshots'
        
        for frame_num in range(self.n_snapshots):
            
            print 'Draw frame number ' + str(frame_num+1)
            snap = self.History[frame_num]        
        
            plt = self.KMC_lat.PlotLattice()            # plot the lattice in this frame
        
            
            for ind in range(self.Species['n_surf']):
                
                # Find all coordinates with species ind occupying it
                x_list = []
                y_list = []
                for site_ind in range(self.n_sites):      # include empty sites
                    if snap[site_ind,2] == ind+1:
                        x_list.append(cart_coords[site_ind,0])
                        y_list.append(cart_coords[site_ind,1])
        
                x = np.array(x_list)
                y = np.array(y_list)                
                
                plt.plot(x, y, linestyle='None', marker = 'o', color = spec_color_list[ind % len(spec_color_list)], markersize = 3, label=spec_label_list[ind])
            
            plt.title('Time: ' + str(self.snap_times[frame_num]) + ' sec')
            plt.legend(frameon=False, loc=4)
                
            plt.savefig(os.path.join(frame_fldr, 'Snapshot_' + str(frame_num+1)))
            plt.close()