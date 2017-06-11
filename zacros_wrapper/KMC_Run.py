import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# For executable
import os
import sys
import subprocess
import copy
import re as _re

from Helper import *
from IO_data import *
from Lattice import Lattice

class kmc_traj():
    
    '''
    Handles a single Zacros trajectory.
    '''
    
    def __init__(self, path = None):
        
        '''
        Initializes class variables
        '''

        self.Path = path
        self.exe_file = None
        
        # Input variables
        self.simin = SimIn()
        self.mechin = MechanismIn()
        self.clusterin = ClusterIn()
        self.lat = Lattice()
        self.statein = StateIn()
        
        # Output variables
        self.genout = PerformanceOut()
        self.specnumout = SpecnumOut()
        self.procstatout = ProcstatOut()
        self.histout = HistoryOut()
        
        # Extra analysis variables
        self.gas_prod = None
        
    '''
    ======================================= File input/output methods =======================================
    '''
    
    def ReadAllInput(self):
    
        '''
        Read all Zacros input files
        '''
        
        self.simin.ReadIn(self.Path)
        self.mechin.ReadIn(self.Path)
        self.clusterin.ReadIn(self.Path)
        self.lat.ReadIn(self.Path)
        self.statein.ReadIn(self.Path)

            
    def WriteAllInput(self):
    
        '''
        Write all Zacros input files
        '''
        
        self.simin.WriteIn(self.Path)
        self.mechin.WriteIn(self.Path)
        self.clusterin.WriteIn(self.Path)
        self.lat.WriteIn(self.Path)
        self.statein.WriteIn(self.Path, self.simin.surf_spec)

        
    def ReadAllOutput(self, build_lattice=False):
        '''
        Read all Zacros output files
        Set build_lattice = True if you want to build the lattice object file.
        This will make it take a lot longer to read
        '''
        self.ReadAllInput()

        if self.CheckComplete():

            # Standard output files
            self.genout.ReadOut(self.Path, self.simin.surf_spec, self.simin.gas_spec )
            self.specnumout.ReadOut(self.Path)
            self.procstatout.ReadOut(self.Path)
            
            if build_lattice:
                self.lat.Read_lattice_output(self.Path)
                
            # Count the number of sites
            with open( os.path.join(self.Path, 'lattice_output.txt'), 'r') as txt:
                RawTxt = txt.readlines()
            nSites = len(RawTxt) - 2
            
            self.histout.ReadOut(self.Path, nSites)
            
            # Extra binary files
            if os.path.isfile(os.path.join(self.Path, 'Prop_output.bin')):
                self.ReadProp(0)
            if os.path.isfile(os.path.join(self.Path, 'PropCounter_output.bin')):
                self.ReadProp(1)
            if os.path.isfile(os.path.join(self.Path, 'SA_output.bin')):
                self.ReadSA()

        else:
            print 'general_output.txt not found in ' + self.Path
    
    
    def ReadProp(self, Mode):
        '''
        Mode = 0: Read Prop_output.bin
        Mode = 1: Read PropCounter_output.bin
        '''
        dt = np.dtype(np.float64)
        if Mode == 0:     # Instantaneous propensities
            FileName = 'Prop_output.bin'
        elif Mode == 1:   # Integral propensities
            FileName = 'PropCounter_output.bin'

        virtual_arr = np.memmap(os.path.join(self.Path, FileName), dt, "r")
        nRxn = len( self.genout.RxnNameList )
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nRxn)
        virtual_arr = virtual_arr[:nNum]

        if Mode == 0:
            self.prop = np.reshape(virtual_arr, [nNum/nRxn, nRxn])
            self.prop = np.array(self.prop[::self.procstatout.Spacing])
            
        if Mode == 1:
            self.propCounter = np.reshape(virtual_arr, [nNum/nRxn, nRxn])
            self.propCounter = np.array(self.propCounter[::self.procstatout.Spacing])
        del virtual_arr

    def ReadSA(self):
        '''
        Read SA_output.bin
        '''
        dt = np.dtype(np.float64)
        FileName = 'SA_output.bin'
        if os.path.isfile(os.path.join(self.Path, FileName)):
            virtual_arr = np.memmap(os.path.join(self.Path,
                                                   FileName), dt, "r")
            nRxn = len(self.genout.Nu)
            nNum = virtual_arr.shape[0]
            nNum = nNum - (nNum % nRxn)
            virtual_arr = virtual_arr[:nNum]

            self.W_sen_anal = np.reshape(virtual_arr, [nNum/nRxn, nRxn])
            self.W_sen_anal = np.array(self.W_sen_anal[::self.specnumout.Spacing])

            del virtual_arr
        else:
            print 'No sensitivity analysis output file'
            
            
    def ReadCluster(self):
        '''
        Read clusterocc.bin
        '''
        dt = np.dtype(np.int32)
        virtual_arr = np.memmap(os.path.join(self.Path, 'clusterocc.bin'), dt, "r")
        nCluster = self.clusterin.get_num_clusters()
        nNum = virtual_arr.shape[0]
        nNum = nNum - (nNum % nCluster)
        self.cluster = np.array(np.reshape(virtual_arr,
                                                    [nNum/nCluster, nCluster])
                                        [::self.specnumout.Spacing])
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
       
    
    def CheckComplete(self):
    
        '''
        Check to see if a Zacros run has completed successfully
        '''
        
        Complete = False
        if os.path.isfile(os.path.join(self.Path, 'general_output.txt')):
            with open(os.path.join(self.Path, 'general_output.txt'), 'r') as txt:
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
        for i in self.mechin.rxn_list:
            for j in i.variant_list:
                j.pre_expon = j.pre_expon * delta_sdf[rxn_ind]
                j.scaledown_factor = j.scaledown_factor * delta_sdf[rxn_ind]
                rxn_ind += 1 
    
    
    def time_search(self, t):
        
        '''
        Given a time, look up the index of the smallest time greater than or equal to that time
        '''
        
        if t > self.specnumout.t[-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind = 0
        while self.specnumout.t[ind] < t:
            ind += 1
            
        return ind
        
        
    def time_search_interp(self, t):
        
        '''
        Get the information necessary to linearly interpolate data between time points
        '''
        
        if t > self.specnumout.t[-1] or t < 0:
            raise Exception('Time is out of range.')
        
        ind_geq = 0
        while self.specnumout.t[ind_geq] < t:
            ind_geq += 1
            
        ind_leq = len( self.specnumout.t ) - 1        # set it initially equal to the last index
        while self.specnumout.t[ind_leq] > t:
            ind_leq -= 1
        
        if ind_geq - ind_leq < 0 or ind_geq - ind_leq > 1:
            raise Exception('Time indices are wrong')
        
        low_frac = 1.0
        if not (ind_geq == ind_leq):
            low_frac = (self.specnumout.t[ind_geq] - t) / (self.specnumout.t[ind_geq] - self.specnumout.t[ind_leq])
            
        high_frac = 1.0 - low_frac
        
        return [[ind_leq, ind_geq], [low_frac, high_frac]]
    
        
    @staticmethod
    def time_sandwich(run1, run2):
        
        '''
        Take two trajectories and append them
        '''
        
        sandwich = copy.deepcopy(run1)
        sandwich.genout.t_final = run1.genout.t_final + run2.genout.t_final
        sandwich.genout.events_occurred = run1.genout.events_occurred + run2.genout.events_occurred
        sandwich.genout.CPU_time = run1.genout.CPU_time + run2.genout.CPU_time
        
        sandwich.specnumout.t = np.concatenate([run1.specnumout.t, run2.specnumout.t[1::] + run1.specnumout.t[-1] * np.ones( len(run2.specnumout.t)-1 )])
        
        n_surf_specs = len( self.simin.surf_spec )
        run2.specnumout.spec[1::, n_surf_specs : ] = run2.specnumout.spec[1::, n_surf_specs : ] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.specnumout.spec[-1, n_surf_specs : ]] )        
        sandwich.specnumout.spec = np.vstack([run1.specnumout.spec, run2.specnumout.spec[1::,:] ])
        sandwich.prop = np.vstack([run1.prop, run2.prop[1::,:] ])
        sandwich.procstatout.events = np.vstack( [run1.procstatout.events, run2.procstatout.events[1::,:] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.procstatout.events[-1,:]] ) ] )
        sandwich.propCounter = np.vstack( [run1.propCounter, run2.propCounter[1::,:] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.propCounter[-1,:]]  ) ] )
        sandwich.W_sen_anal  = np.vstack( [run1.W_sen_anal, run2.W_sen_anal[1::,:] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.W_sen_anal[-1,:]]  ) ] )      
        
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
        for i in range (len( self.simin.surf_spec )):
            time_vecs.append(self.specnumout.t)
            surf_spec_vecs.append(self.specnumout.spec[:,i] / float(site_norm))
        
        PlotTimeSeries(time_vecs, surf_spec_vecs, xlab = 'Time (s)', ylab = ylabel, series_labels = self.simin.surf_spec, fname = os.path.join(self.Path, 'surf_spec_vs_time.png'))
    
    
    def PlotGasSpecVsTime(self):
        
        '''
        Plot gas phase species profiles versus time
        '''
        
        time_vecs = []
        gas_spec_vecs = []
        for i in range(len( self.simin.gas_spec )):
            time_vecs.append(self.specnumout.t)
            gas_spec_vecs.append(self.specnumout.spec[:, i + len( self.simin.surf_spec ) ])
        
        PlotTimeSeries(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'Spec. pop.', series_labels = self.simin.gas_spec, fname = os.path.join(self.Path, 'gas_spec_vs_time.png'))
        
        
    
    def PlotElemStepFreqs(self, window = [0.0, 1.0], time_norm = False, site_norm = 1):
        
        '''
        Plot a bar graph of elementary step frequencies versus time
        '''
        
        start_ind = self.time_search(window[0] * self.specnumout.t[-1])
        end_ind = self.time_search(window[1] * self.specnumout.t[-1])
        event_freqs = ( self.procstatout.events[end_ind,:] - self.procstatout.events[start_ind,:] ) / float(site_norm)
        if time_norm:
            event_freqs = event_freqs / ( self.specnumout.t[end_ind] - self.specnumout.t[start_ind] )
        
        PlotOptions()
        plt.figure()
        
        width = 0.2
        ind = 0
        yvals = []
        ylabels = []
        bar_vals = []
        
        store_ind = 0       # index of where the data is stored
        
        for rxn in self.mechin.rxn_list:
            for varnt in rxn.variant_list:
        
                fwd_rate = event_freqs[store_ind]
                store_ind += 1
                bwd_rate = 0
                
                if rxn.is_reversible:
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
                    ylabels.append( rxn.name + varnt.name )
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
        plt.tight_layout()
        
        plt.savefig(os.path.join(self.Path, 'elem_step_freqs.png'))
        plt.close()
        
        
    def PlotLattice(self):
    
        '''
        Plot the lattice
        '''
        
        if self.lat.text_only:
        
            print 'Cannot plot lattice. Only text input exists.'
            
        else:
        
            plt = self.lat.PlotLattice()
            plt.savefig(os.path.join(self.Path, 'lattice.png'))
            plt.close()
        
        
    def LatticeMovie(self, include_neighbor_lines = False, spec_color_list = ['b', 'g','r','c','m','y','k']):       # Need make marker type consistent with the site type

        '''
        Create a .png file with a picture of a lattice for every snapshot in history_output.txt
        '''
    
        cart_coords = self.KMC_lat.cart_coords
        spec_label_list = self.simin.surf_spec
        
        frame_fldr = os.path.join(self.Path, 'lattice_frames')
        if not os.path.exists( frame_fldr ):
                os.makedirs( frame_fldr )
                
        print str(self.histout.n_snapshots) + ' total snapshots'
        
        for frame_num in range(self.histout.n_snapshots):
            
            print 'Draw frame number ' + str(frame_num+1)
            snap = self.histout.snapshots[frame_num]
        
            plt = self.KMC_lat.PlotLattice()            # plot the lattice in this frame
        
            
            for ind in range( len( self.simin.surf_spec ) ):
                
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