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

class kmc_traj():
    
    '''
    Handles a single Zacros trajectory.
    '''
    
    def __init__(self):
        
        '''
        Initializes class variables
        '''

        self.exe_file = None
        self.props_avg = None
        self.rate_traj = None
        self.int_rate_traj = None
        self.gas_prod = None
        
    '''
    ======================================= File input/output methods =======================================
    '''
    
    def ReadAllInput(self):
        '''
        Read all input files. state_input.dat will be read only
        if it is present.
        '''
        self.ReadSimIn()
        self.ReadLatticeIn()
        self.ReadEngIn()
        self.ReadMechIn()

        if _os.path.isfile(_os.path.join(self.Path, 'state_input.dat')):
            self.ReadStateInput()
        else:
            self.State = None

            
    def WriteAllInput(self):
        '''
        Write all Zacros input files
        '''
        self.WriteSimIn()
        self.WriteMechanism()
        self.WriteEnergetics()
        self.WriteStateIn()
        self.WriteLattice()

        
    def ReadAllOutput(self, build_lattice=False):
        '''
        Read all Zacros output files
        Set build_lattice = True if you want to build the lattice object file.
        This will make it take a lot longer to read
        '''
        self.ReadAllInput()

        if self.CheckComplete():

            # Standard output files
            self.ReadGeneral()
            self.ReadProcstat()
            self.ReadSpecnum()
            if build_lattice:
                self.KMC_lat.Read_lattice_output(_os.path.join
                                                 (self.Path,
                                                  'lattice_output.txt'))
            self.ReadHistory()

            # Extra binary files
            if _os.path.isfile(_os.path.join(self.Path,
                                             'Prop_output.bin')):
                self.ReadProp(0)
            if _os.path.isfile(_os.path.join(self.Path,
                                             'PropCounter_output.bin')):
                self.ReadProp(1)
            if _os.path.isfile(_os.path.join(self.Path,
                                             'SA_output.bin')):
                self.ReadSA()

        else:
            print 'general_output.txt not found in ' + self.Path
    
    
    def ReadProp(self, Mode):
        '''
        Mode = 0: Read Prop_output.bin
        Mode = 1: Read PropCounter_output.bin
        '''
        dt = _np.dtype(_np.float64)
        if Mode == 0:     # Instantaneous propensities
            FileName = 'Prop_output.bin'
        elif Mode == 1:   # Integral propensities
            FileName = 'PropCounter_output.bin'

        virtual_arr = _np.memmap(_os.path.join(self.Path,
                                               FileName), dt, "r")
        nRxn = len(self.Performance.Nu)
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]
        if not hasattr(self, 'Binary'):
            self.Binary = BinaryIn()

        if Mode == 0:
            self.Binary.prop = _np.reshape(virtual_arr, [nNum/nRxn, nRxn])
            self.Binary.prop = _np.array(self.Binary.prop
                                         [::self.Procstat.Spacing])
        if Mode == 1:
            self.Binary.propCounter = _np.reshape(virtual_arr,
                                                  [nNum/nRxn, nRxn])
            self.Binary.propCounter = _np.array(self.Binary.propCounter
                                                [::self.Procstat.Spacing])
        del virtual_arr

    def ReadSA(self):
        '''
        Read SA_output.bin
        '''
        dt = _np.dtype(_np.float64)
        FileName = 'SA_output.bin'
        if _os.path.isfile(_os.path.join(self.Path, FileName)):
            virtual_arr = _np.memmap(_os.path.join(self.Path,
                                                   FileName), dt, "r")
            nRxn = len(self.Performance.Nu)
            nNum = virtual_arr.shape[0]
            nNum = nNum - (nNum % nRxn)
            virtual_arr = virtual_arr[:nNum]
            if not hasattr(self, 'Binary'):
                self.Binary = BinaryIn()

            self.Binary.W_sen_anal = _np.reshape(virtual_arr,
                                                 [nNum/nRxn, nRxn])
            self.Binary.W_sen_anal = _np.array(self.Binary.W_sen_anal
                                               [::self.Specnum.Spacing])

            del virtual_arr
        else:
            print 'No sensitivity analysis output file'
            
            
    def ReadCluster(self):
        '''
        Read clusterocc.bin
        '''
        dt = _np.dtype(_np.int32)
        virtual_arr = _np.memmap(_os.path.join(self.Path, 'clusterocc.bin'),
                                 dt, "r")
        nCluster = sum(len(s.variant_name) for s in self.Cluster)
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nCluster)
        if not hasattr(self, 'Binary'):
            self.Binary = BinaryIn()
        self.Binary = BinaryIn()
        self.Binary.cluster = _np.array(_np.reshape(virtual_arr,
                                                    [nNum/nCluster, nCluster])
                                        [::self.Specnum.Spacing])
        del virtual_arr
    
    '''
    ======================================= Calculation methods =======================================
    '''
    
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
       
    
    def CheckComplete(self):        # Copy/pasted from IOdata, needs to be made compatible
        '''
        Check to see if a Zacros run has completed successfully
        '''
        Complete = False
        if _os.path.isfile(_os.path.join(self.Path,
                                         'general_output.txt')):
            with open(_os.path.join(self.Path,
                                    'general_output.txt'),
                      'r') as txt:
                RawTxt = txt.readlines()
            for i in RawTxt:
                if _re.search('Normal termination', i):
                    Complete = True
        return Complete
        
    
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
    
    
    def time_search(self, t):
        
        '''
        Given a time, look up the index of the smallest time greater than or equal to that time
        '''
        
        if t > self.Specnum['t'][-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind = 0
        while self.Specnum['t'][ind] < t:
            ind += 1
            
        return ind
        
        
    def time_search_interp(self, t):
        
        '''
        Get the information necessary to linearly interpolate data between time points
        '''
        
        if t > self.Specnum['t'][-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind_geq = 0
        while self.Specnum['t'][ind_geq] < t:
            ind_geq += 1
            
        ind_leq = len(self.Specnum['t']) - 1        # set it initially equal to the last index
        while self.Specnum['t'][ind_leq] > t:
            ind_leq -= 1
        
        if ind_geq - ind_leq < 0 or ind_geq - ind_leq > 1:
            raise Exception('Time indices are wrong')
        
        low_frac = 1.0
        if not (ind_geq == ind_leq):
            low_frac = (self.Specnum['t'][ind_geq] - t) / (self.Specnum['t'][ind_geq] - self.Specnum['t'][ind_leq])
            
        high_frac = 1.0 - low_frac
        
        return [[ind_leq, ind_geq], [low_frac, high_frac]]
    
        
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
        sandwich.Binary['prop'] = np.vstack([run1.Binary['prop'], run2.Binary['prop'][1::,:] ])
        sandwich.Procstat['events'] = np.vstack( [run1.Procstat['events'], run2.Procstat['events'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Procstat['events'][-1,:]] ) ] )
        sandwich.Binary['propCounter'] = np.vstack( [run1.Binary['propCounter'], run2.Binary['propCounter'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Binary['propCounter'][-1,:]]  ) ] )
        sandwich.Binary['W_sen_anal']  = np.vstack( [run1.Binary['W_sen_anal'], run2.Binary['W_sen_anal'][1::,:] + np.dot(np.ones([len(run2.Specnum['t'])-1 ,1]), [run1.Binary['W_sen_anal'][-1,:]]  ) ] )      
        
        sandwich.History = [ run1.History[0], run2.History[-1] ]     
        
        return sandwich
    

    '''
    ==================================== Plotting methods ====================================
    '''
    
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
        
        start_ind = self.time_search(window[0] * self.Specnum['t'][-1])
        end_ind = self.time_search(window[1] * self.Specnum['t'][-1])
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
        
        
    def LatticeMovie(self, include_neighbor_lines = False, spec_color_list = ['b', 'g','r','c','m','y','k']):       # Need make marker type consistent with the site type

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