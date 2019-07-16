import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches

# For executable
import os
import sys
import subprocess
import copy
import re as _re

from zacros_wrapper.utils import *
from zacros_wrapper.IO_data import *
from zacros_wrapper.Lattice import Lattice

class kmc_traj():

    '''
    Handles a single Zacros trajectory.
    '''

    def __init__(self, path = None):

        '''
        Initialize class variables
        '''

        self.Path = path                        # Path to the folder with the Zacros input and output files
        self.exe_file = None                    # Path to the Zacros executable

        # Input variables
        self.simin = SimIn()                    # data from simulation_input.dat
        self.mechin = MechanismIn()             # data from mechanism_input.dat
        self.clusterin = ClusterIn()            # data from energetics_input.dat
        self.lat = Lattice()                    # data from lattice_input.dat
        self.statein = StateIn()                # data from state_input.dat

        # Output variables
        self.genout = PerformanceOut()          # data from general_output.txt
        self.specnumout = SpecnumOut()          # data from specnum_output.txt
        self.procstatout = ProcstatOut()        # data from procstat_output.txt
        self.histout = HistoryOut()             # data from history_output.txt

        # Extra analysis variables
        self.gas_prod = None                    # string with the name of the product (gas) species, used to compute the rate

        # Missing other variables obtained from additional output files...??

        # Additional output
        self.prop = None                        # data from propensities_output.txt
        self.propCounter = None                 # data from timeintprop_output.txt
        self.W_sen_anal = None                  # data from trajderiv_output.txt
        self.spec_num_int = None                # data from timeintspecs_output.txt
        #self.TS_site_props_list = None
        #self.TS_site_props_ss = None            # steady state site propensities


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

        :param build_lattice :     True - builds a Lattice object
            False - reads lattice_output.txt as text only
        '''

        #print 'Reading kMC trajectory data from ' + self.Path

        self.ReadAllInput()

        if self.CheckComplete():

            # Standard output files
            self.genout.ReadOut(self.Path, self.simin.surf_spec, self.simin.gas_spec )
            self.specnumout.ReadOut(self.Path)
            if os.path.isfile(os.path.join(self.Path, 'procstat_output.txt')):
                self.procstatout.ReadOut(self.Path)

            if build_lattice:
                self.lat.Read_lattice_output(self.Path)

            # Count the number of sites
            with open( os.path.join(self.Path, 'lattice_output.txt'), 'r') as txt:
                RawTxt = txt.readlines()
            nSites = len(RawTxt) - 2

            self.histout.ReadOut(self.Path, nSites)

            # Extra output files
			# Disable for simplicity
            #self.prop = Read_propensities(self.Path, len( self.genout.RxnNameList ) )
            #self.propCounter = Read_time_integrated_propensities(self.Path, len( self.genout.RxnNameList ) )
            #self.W_sen_anal = Read_trajectory_derivatives(self.Path, len( self.genout.RxnNameList ))
            #self.spec_num_int = Read_time_integrated_species(self.Path, len( self.simin.surf_spec ))
            #self.TS_site_props_list = Read_time_integrated_site_props(self.Path, nSites, len( self.genout.RxnNameList ), self.histout.n_snapshots )


            #if not self.TS_site_props_list is None:

                #if self.histout.snap_times[-1] == 0.:   # only 1 entry in history output
                    #self.TS_site_props_ss = self.TS_site_props_list[-1]
                #else:
                    #self.TS_site_props_ss = ( self.TS_site_props_list[-1] - self.TS_site_props_list[0] ) / ( self.histout.snap_times[-1] - self.histout.snap_times[0] )

        else:
            print('general_output.txt not found in ' + self.Path)


    '''
    ======================================= Calculation methods =======================================
    '''


    def Run_sim(self):

        '''
        Call the Zacros executable and run the simulation. self.exe_file must have been initialized
        with the full file path of the Zacros executable.
        '''

        os.chdir(self.Path)

        if self.exe_file is None:
            raise Exception('Zacros executable not specified.')

        try:
            print('--- Zacros run starting ---')
            subprocess.call([self.exe_file])
            print('--- Zacros run completed ---')
        except:
            raise Exception('Zacros run in ' + self.Path + ' failed.')


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

        :param delta_sdf: list or vector of ratios to apply. Must have length equal to the number of reactions.
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

        :param t:     time to search for
        '''

        if t > self.specnumout.t[-1] or t < 0:
            raise Exception('Time is out of range.')

        ind = 0
        while self.specnumout.t[ind] < t:
            ind += 1

        return ind


    def time_avg_covs(self, t1 = 0, t2 = None):
        '''
        Time average surface species counts

        :param t1: Start time for time averaging interval
        :param t2: End time for time averaging interval, if left as None it will be set as the final time
        '''

        if t2 is None:          # Use final time by default
            t2 = self.specnumout.t[-1]

        idb_start = self.time_search_interp(t1)
        idb_end = self.time_search_interp(t2)

        cov_integ_start = idb_start[1][0] * self.spec_num_int[idb_start[0][0], :] + idb_start[1][1] * self.spec_num_int[idb_start[0][1], :]
        cov_integ_end = idb_end[1][0] * self.spec_num_int[idb_end[0][0], :] + idb_end[1][1] * self.spec_num_int[idb_end[0][1], :]

        return ( cov_integ_end  - cov_integ_start ) / (t2 - t1)

    def time_search_interp(self, t):

        '''
        Get the information necessary to linearly interpolate data between time points

        :param t: Time at which you want to compute interpolated data
        :returns: A 2-item list of 2-items each in the order: Index of time point less than or equal to t,
            index of time point greater than or equal to t
            weighting factor for lower point
            weighting factor for higher point
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


    '''
    ==================================== Plotting methods ====================================
    '''

    def PlotSurfSpecVsTime(self, site_norm = 1, savefig = True, incl_total = False):

        '''
        Plot surface species profiles versus time - output in surf_spec_vs_time.png in the directory with the Zacros run

        :param site_norm: Normalizing factor for the coverages (e.g. nubmer of top sites). Default value if 1.
        '''

        if site_norm == 1:
            ylabel = 'Species count'
        else:
            ylabel = 'Coverage'

        time_vecs = []
        surf_spec_vecs = []
        surf_spec_name = []

        for i in range (len( self.simin.surf_spec )):

            if np.any(self.specnumout.spec[:,i]): # Pick all the non-zero (observable) speices

                time_vecs.append(self.specnumout.t)
                surf_spec_vecs.append(self.specnumout.spec[:,i] / float(site_norm))
                surf_spec_name.append(self.simin.surf_spec[i])

        if incl_total:

            time_vecs.append(self.specnumout.t)
            surf_spec_vecs.append(np.sum(surf_spec_vecs, axis = 0))
            surf_spec_name.append('Total')

        if savefig:

            PlotTimeSeries(time_vecs, surf_spec_vecs, xlab = 'Time (s)', ylab = ylabel, series_labels = surf_spec_name, fname = os.path.join(self.Path, 'surf_spec_vs_time.png'))

        else:

            PlotTimeSeries(time_vecs, surf_spec_vecs, xlab = 'Time (s)', ylab = ylabel, series_labels = surf_spec_name)


    def PlotGasSpecVsTime(self, savefig = True):

        '''
        Plot gas phase species profiles versus time - output in gas_spec_vs_time.png in the directory with the Zacros run
        '''

        time_vecs = []
        gas_spec_vecs = []
        for i in range(len( self.simin.gas_spec )):
            time_vecs.append(self.specnumout.t)
            gas_spec_vecs.append(self.specnumout.spec[:, i + len( self.simin.surf_spec ) ])

        if savefig:

            PlotTimeSeries(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'Spec. pop.', series_labels = self.simin.gas_spec, fname = os.path.join(self.Path, 'gas_spec_vs_time.png'))

        else:

            PlotTimeSeries(time_vecs, gas_spec_vecs, xlab = 'Time (s)', ylab = 'Spec. pop.', series_labels = self.simin.gas_spec)


    def PlotElemStepFreqs(self, window = [0.0, 1.0], time_norm = False, site_norm = 1, savefig = True):

        '''
        Plot a bar graph of elementary step frequencies versus time - output in elem_step_freqs.png in the directory with the Zacros run

        :param window: Beginning and ending fraction of the trajectory over which to count steps.

        :param time_norm: True - normalizes event frequencies by time. False - Plots total event firings.

        :param site_norm: Normalization factor for event frequencies.
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
        plt.yticks(size=12)
        plt.xlabel('Frequency',size=24)
        plt.yticks(yvals, ylabels)
        plt.xlim([xmin, xmax])


        plt.legend(['fwd', 'rev', 'net'], bbox_to_anchor = (1.05, 1),loc= 'upper left', prop={'size':12},frameon=False)

        ax = plt.gca()
        leg = ax.get_legend()
        leg.legendHandles[0].set_color('red')
        leg.legendHandles[1].set_color('blue')
        leg.legendHandles[2].set_color('green')

        plt.tight_layout()

        if savefig:
            plt.savefig(os.path.join(self.Path, 'elem_step_freqs.png'), bbox_inches = "tight")

        else:
            plt.show()

        plt.close()


    def PlotLattice(self, savefig = True):

        '''
        Plot the lattice - output in lattice.png in the run directory
        '''

        if self.lat.text_only:

            print('Cannot plot lattice. Only text input exists.')

        else:
            plt = self.lat.PlotLattice()

            if savefig:
                plt.savefig(os.path.join(self.Path, 'lattice.png'))
            else:
                plt.show()
            plt.close()

    def PlotLattice(self, savefig = True):

        '''
        Plot the lattice - output in lattice.png in the run directory
        '''

        if self.lat.text_only:

            print('Cannot plot lattice. Only text input exists.')

        else:
            plt = self.lat.PlotLattice()

            if savefig:
                plt.savefig(os.path.join(self.Path, 'lattice.png'))
            else:
                plt.show()
            plt.close()

    def PlotLattice3D(self, dz = 1, get_GIF = False, savefig = True):
        '''
        Plot the 3D lattice - output in - output in lattice.png in the run directory
        '''

        if self.lat.text_only:

            print('Cannot plot lattice. Only text input exists.')

        else:
            self.lat.set_cart_coords_3D(dz)
            fig, _ = self.lat.PlotLattice3D(get_GIF = get_GIF, type_symbols=['.'], selected_sites = 'l0')

            if not get_GIF:
                if savefig:
                    fig.savefig(os.path.join(self.Path, 'lattice.png'))
                else:
                    fig.show()
                plt.close()

    def LatticeMovie3D(self, dz = 1, include_neighbor_lines = False, spec_color_list = ['slateblue', 'green', 'maroon','peru'], savefig = True):       # Need make marker type consistent with the site type

        '''
        Create a subfolder called lattice_frames
        Create a .png file with a picture of a lattice for every snapshot in history_output.txt

        include_neighbor_lines :    If true, will draw lines between neighboring lattice sites (takes more time)
        spec_color_list :           List of colors to use for different species, will cycle through if there are more species than colors
        '''
        self.lat.set_cart_coords_3D(dz)
        cart_coords_3d = self.lat.cart_coords_3d
        spec_label_list = self.simin.surf_spec

        frame_fldr = os.path.join(self.Path, 'lattice_frames_3D')
        if not os.path.exists( frame_fldr ):
                os.makedirs( frame_fldr )

        print(str(self.histout.n_snapshots) + ' total snapshots')

        for frame_num in range(int(self.histout.n_snapshots)):

            print('Draw frame number ' + str(frame_num+1))
            snap = self.histout.snapshots[frame_num]

            fig, ax = self.lat.PlotLattice3D(type_symbols=['.'],  selected_sites = 'l0')            # plot the lattice in this frame

            for ind in range( len( self.simin.surf_spec ) ):

                # Find all coordinates with species ind occupying it
                x_list = []
                y_list = []
                z_list = []

                for site_ind in range(cart_coords_3d.shape[0]):      # include empty sites
                    if snap[site_ind,2] == ind+1:
                        x_list.append(cart_coords_3d[site_ind,0])
                        y_list.append(cart_coords_3d[site_ind,1])
                        z_list.append(cart_coords_3d[site_ind,2])

                x = np.array(x_list)
                y = np.array(y_list)
                z = np.array(z_list)

                for xi, yi, zi in zip(x,y,z):
                    layer_i = int(zi -1)
                    ax.scatter3D(xi, yi, zi, marker = 'o', color = spec_color_list[int(layer_i% len(spec_color_list))],  s = 150, edgecolors = 'k')

            # Create legend labels
            n_layers = int(np.max(z))
            legend_patches = []
            for li in range(n_layers):
                layer_name = self.lat.site_type_names[li]
                legend_patches.append(mpatches.Patch(color=spec_color_list[int(li% len(spec_color_list))], label= spec_label_list[ind] + '_' + layer_name))

            plt.legend(handles=legend_patches, bbox_to_anchor = (1.02,1), loc = 'upper left', frameon = False)
            plt.title('Time: ' + str(self.histout.snap_times[frame_num]) + ' sec')

            if savefig:
                fig.savefig(os.path.join(frame_fldr, 'Snapshot_' + str(frame_num+1)), bbox_inches = "tight")
            else:
                fig.show()
            plt.close()

    def PlotInitialState(self, dz = 1, include_neighbor_lines = False, spec_color_list = ['slateblue', 'green', 'maroon','peru'], savefig = True):       # Need make marker type consistent with the site type

        '''
        Create a .png file with a picture of a lattice for the initial snapshot given state_input.dat

        include_neighbor_lines :    If true, will draw lines between neighboring lattice sites (takes more time)
        spec_color_list :           List of colors to use for different species, will cycle through if there are more species than colors
        '''
        # Count the number of sites
        with open( os.path.join(self.Path, 'lattice_output.txt'), 'r') as txt:
            RawTxt = txt.readlines()
        nSites = len(RawTxt) - 2

        # Read state_input.data
        self.statein.ReadSeed(nSites)

        self.lat.set_cart_coords_3D(dz)
        cart_coords_3d = self.lat.cart_coords_3d
        spec_label_list = self.simin.surf_spec

        frame_fldr = self.Path

        # Plot the initial snap shot
        snap = self.statein.initialsnap

        fig, ax = self.lat.PlotLattice3D(type_symbols=['.'],  selected_sites = 'l0')            # plot the lattice in this frame

        for ind in range( len( self.simin.surf_spec ) ):

            # Find all coordinates with species ind occupying it
            x_list = []
            y_list = []
            z_list = []

            for site_ind in range(cart_coords_3d.shape[0]):      # include empty sites
                if snap[site_ind,2] == ind+1:
                    x_list.append(cart_coords_3d[site_ind,0])
                    y_list.append(cart_coords_3d[site_ind,1])
                    z_list.append(cart_coords_3d[site_ind,2])

            x = np.array(x_list)
            y = np.array(y_list)
            z = np.array(z_list)

            for xi, yi, zi in zip(x,y,z):
                layer_i = int(zi -1)
                ax.scatter3D(xi, yi, zi, marker = 'o', color = spec_color_list[int(layer_i% len(spec_color_list))],  s = 150, edgecolors = 'k')

        # Create legend labels
        n_layers = int(np.max(z))
        legend_patches = []
        for li in range(n_layers):
            layer_name = self.lat.site_type_names[li]
            legend_patches.append(mpatches.Patch(color=spec_color_list[int(li% len(spec_color_list))], label= spec_label_list[ind] + '_' + layer_name))

        plt.legend(handles=legend_patches, bbox_to_anchor = (1.02,1), loc = 'upper left', frameon = False)
        plt.title('Time = 0 ' + ' sec')

        if savefig:
            fig.savefig(os.path.join(frame_fldr, 'Initial_Snapshot'), bbox_inches = "tight")
        else:
            fig.show()
        plt.close()


    def LatticeMovie(self, include_neighbor_lines = False, spec_color_list = colors_pool, savefig = True):       # Need make marker type consistent with the site type

        '''
        Create a subfolder called lattice_frames
        Create a .png file with a picture of a lattice for every snapshot in history_output.txt

        include_neighbor_lines :    If true, will draw lines between neighboring lattice sites (takes more time)
        spec_color_list :           List of colors to use for different species, will cycle through if there are more species than colors
        '''

        cart_coords = self.lat.cart_coords
        spec_label_list = self.simin.surf_spec

        frame_fldr = os.path.join(self.Path, 'lattice_frames')
        if not os.path.exists( frame_fldr ):
                os.makedirs( frame_fldr )

        print(str(self.histout.n_snapshots) + ' total snapshots')

        for frame_num in range(int(self.histout.n_snapshots)):

            print('Draw frame number ' + str(frame_num+1))
            snap = self.histout.snapshots[frame_num]

            plt = self.lat.PlotLattice()            # plot the lattice in this frame


            for ind in range( len( self.simin.surf_spec ) ):

                # Find all coordinates with species ind occupying it
                x_list = []
                y_list = []
                for site_ind in range(cart_coords.shape[0]):      # include empty sites
                    if snap[site_ind,2] == ind+1:
                        x_list.append(cart_coords[site_ind,0])
                        y_list.append(cart_coords[site_ind,1])

                x = np.array(x_list)
                y = np.array(y_list)

                plt.plot(x, y, linestyle='None', marker = 'o', color = spec_color_list[ind % len(spec_color_list)], markersize = 5.5, markeredgewidth = 0.0, label=spec_label_list[ind])

            plt.title('Time: ' + str(self.histout.snap_times[frame_num]) + ' sec')
            plt.legend(bbox_to_anchor = (1.02,1), loc = 'upper left', prop = {'size':12}, frameon = False)

            if savefig:
                plt.savefig(os.path.join(frame_fldr, 'Snapshot_' + str(frame_num+1)), bbox_inches = "tight")
            else:
                plt.show()
            plt.close()


def append_trajectories(run1, run2):

    '''
    Take two trajectories and append them together.

    :param run1: First trajectory - a kmc_traj object
    :param run2: Second trajectory - a kmc_traj object which continues the simulation where the first one terminated
    :returns: A kmc_traj object that has the trajectories appended
    '''

    combo = copy.deepcopy(run1)
    combo.genout.t_final = run1.genout.t_final + run2.genout.t_final
    combo.genout.events_occurred = run1.genout.events_occurred + run2.genout.events_occurred
    combo.genout.CPU_time = run1.genout.CPU_time + run2.genout.CPU_time

    combo.specnumout.t = np.concatenate([run1.specnumout.t, run2.specnumout.t[1::] + run1.specnumout.t[-1] * np.ones( len(run2.specnumout.t)-1 )])

    n_surf_specs = len( run1.simin.surf_spec )
    run2.specnumout.spec[1::, n_surf_specs : ] = run2.specnumout.spec[1::, n_surf_specs : ] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.specnumout.spec[-1, n_surf_specs : ]] )
    combo.specnumout.spec = np.vstack([run1.specnumout.spec, run2.specnumout.spec[1::,:] ])
    combo.prop = np.vstack([run1.prop, run2.prop[1::,:] ])
    combo.procstatout.events = np.vstack( [run1.procstatout.events, run2.procstatout.events[1::,:] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.procstatout.events[-1,:]] ) ] )
    combo.propCounter = np.vstack( [run1.propCounter, run2.propCounter[1::,:] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.propCounter[-1,:]]  ) ] )
    combo.W_sen_anal  = np.vstack( [run1.W_sen_anal, run2.W_sen_anal[1::,:] + np.dot(np.ones([len(run2.specnumout.t)-1 ,1]), [run1.W_sen_anal[-1,:]]  ) ] )

    hist_t_total = run1.histout.snap_times[-1] + run2.histout.snap_times[-1]
    #combo.TS_site_props_ss = ( run1.TS_site_props_ss * run1.histout.snap_times[-1] + run2.TS_site_props_ss * run2.histout.snap_times[-1] ) / hist_t_total

    #combo.History = [ run1.History[0], run2.History[-1] ]

    return combo
