import os
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
import copy
import sys
import scipy.stats

from zacros_wrapper.KMC_Run import *
from zacros_wrapper.utils import *
import time
import scipy

class Replicates:

    '''
    Performs statistical data analysis for muliple kMC trajectories with the same input, but different
    random seeds and possibly different initial states
    '''

    def __init__(self):

        '''
        Initialize data structures
        '''

        # General info
        self.ParentFolder = None
        self.runtemplate = None                             # Use this to build replicate jobs
        self.runAvg = None            # Values are averages of all runs from runList
        self.avg_updated = False        # Whether we need to update the trajectory averages
        self.n_trajectories = 0

        # Input data for different trajectories
        self.run_dirs = []
        self.rand_seeds = []

        # Output data taken accross all trajectories
        # Put here various multidimensional numpy arrays which will make it easier to do statistics
        self.t_vec = None                    # times at which data is recorded
        self.species_pops = None            # species populations
        self.rxn_freqs = None               # reaction frequencies
        self.History_final_snaps = None
        self.propensities = None            # total propensity for each reaction
        self.Props_integ = None
        self.traj_derivs = None
        self.events_total = None
        self.CPU_total = None
        self.time_avg_covs = None
        #self.TS_site_props_sss = None

        # Rate analysis
        self.N_batches = 1000       # Used to be 1000
        self.Nbpt = None
        self.batch_length = None    # Length of time in one batch
        self.gas_product = None
        self.gas_prod_ind = None
        self.TOF_stoich = None

        # Rate and autocorrelation data
        self.rate = None
        self.rate_CI = None
        self.ACF = None
        self.ACF_CI = None

        # Sensitivity indices
        self.NSC = None
        self.NSC_ci = None


    def Set_kmc_template(self, kmc_temp):

        '''
        Set the KMC trajectory template as well as other useful variables
        '''

        self.runtemplate = copy.deepcopy(kmc_temp)
        self.gas_product = kmc_temp.gas_prod


    def BuildJobFiles(self, init_states = None):

        '''
        Builds folders with Zacros input files. Each trajectory is assigned a different random seed.
        :param init_states: List of intial states for each trajectory.
        '''

        if not os.path.exists(self.ParentFolder):
            os.makedirs(self.ParentFolder)
        ClearFolderContents(self.ParentFolder)

        # List the directories and random seeds
        self.run_dirs = []
        self.rand_seeds = []
        seed = 5000

        # Go through the directories and write input files for trajectories with
        # different random seeds and possibly different initial states
        for i in range(self.n_trajectories):

            # Set the random seed
            self.runtemplate.simin.Seed = seed
            self.rand_seeds.append(seed)
            seed = seed + 1

            # Set the path
            fldr = os.path.join(self.ParentFolder, str(i+1))
            self.runtemplate.Path = fldr
            self.run_dirs.append(fldr)

            # Set initial state
            if not init_states is None:
                self.runtemplate.statein.Type = 'history'
                self.runtemplate.statein.Struct = init_states[i]

            # Write input files in that folder
            if not os.path.exists(fldr):
                os.makedirs(fldr)

            self.runtemplate.WriteAllInput()


    def RunAllTrajectories_JobArray(self, max_cores = 100, server = 'Squidward', job_name = 'zacros_JA'):

        '''
        Runs a job array on Squidward or Farber
        Writes the .qs file with certain parameters changed
        Then, it submits it to the gridengine and waits for them to finish
        :param max_cores: Maximum number of cores to use
        :param server: Name of the server. Squidward and Farber are supported.
        :param job_name: Name of the job to put in the submit file.
        '''

        sys.stdout.write('Running parallel jobs\n')
        sys.stdout.flush()

        with open(os.path.join(self.ParentFolder, 'dir_list.txt'), 'w') as txt:
            for fldr in self.run_dirs:
                txt.write(fldr + '\n')

        if server == 'Squidward':
            max_cores = np.min([max_cores, 96])
        elif server == 'Farber':
            max_cores = np.min([max_cores, 100])
        else:
            raise Exception('Unrecognized server for parallel runs')

        n_cores = np.min([max_cores, self.n_trajectories])

        with open(os.path.join(self.ParentFolder, 'zacros_submit_JA.qs'), 'w') as txt:

            if server == 'Farber':

                txt.write('#!/bin/bash\n')
                txt.write('#$ -cwd\n')
                txt.write('#$ -j y\n')
                txt.write('#$ -S /bin/bash\n')
                txt.write('#$ -l h_cpu=168:00:00\n')
                txt.write('#\n')
                txt.write('\n')
                txt.write('#$ -N ' + job_name + ' 					#This is the name of the job array\n')
                txt.write('#$ -t 1-' + str(self.n_trajectories) + '  							#Assumes task IDs increment by 1; can also increment by another value\n')
                txt.write('#$ -tc ' + str(n_cores) + ' 							#This is the total number of tasks to run at any given moment\n')
                txt.write('#$ -pe threads 1 -l mem_free=2G              #Change the last field to the number of processors desired per task\n')
                txt.write('#\n')
                txt.write('# Change the following to #$ and set the amount of memory you need\n')
                txt.write('# per-slot if you are getting out-of-memory errors using the\n')
                txt.write('# default:\n')
                txt.write('#$ -l m_mem_free=4G\n')
                txt.write('\n')
                txt.write('source /etc/profile.d/valet.sh\n')
                txt.write('\n')
                txt.write('# Use vpkg_require to setup the environment:\n')
                txt.write('vpkg_require intel/2016\n')
                txt.write('\n')
                txt.write('# Ensure that the OpenMP runtime knows how many processors to use;\n')
                txt.write('# Grid Engine automatically sets NSLOTS to the number of cores granted\n')
                txt.write('# to this job:\n')
                txt.write('export OMP_NUM_THREADS=$NSLOTS\n')
                txt.write('\n')
                txt.write('job_file=\'' + os.path.join(self.ParentFolder, 'dir_list.txt') + '\'\n')
                txt.write('#Change to the job directory\n')
                txt.write('job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")\n')
                txt.write('cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>\n')
                txt.write('                  #This could easily be modified to take a prefix; ask me how.\n')
                txt.write('\n')
                txt.write('# Now append whatever commands you use to run your OpenMP code:\n')
                txt.write('time ' + self.runtemplate.exe_file)

            elif server == 'Squidward':

                txt.write('#!/bin/bash\n')
                txt.write('#$ -cwd\n')
                txt.write('#$ -j y\n')
                txt.write('#$ -S /bin/bash\n')
                txt.write('#\n')
                txt.write('\n')
                txt.write('#$ -N ' + job_name + ' 					#This is the name of the job array\n')
                txt.write('#$ -t 1-' + str(self.n_trajectories) + '  							#Assumes task IDs increment by 1; can also increment by another value\n')
                txt.write('#$ -tc ' + str(n_cores) + ' 							#This is the total number of tasks to run at any given moment\n')
                txt.write('#$ -pe openmpi-smp 1 -l mem_free=1G			#Change the last field to the number of processors desired per task\n')
                txt.write('\n')
                txt.write('job_file=\'' + os.path.join(self.ParentFolder, 'dir_list.txt') + '\'\n')
                txt.write('#Change to the job directory\n')
                txt.write('job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")\n')
                txt.write('cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>\n')
                txt.write('\n\n')
                txt.write('time ' + self.runtemplate.exe_file)

            else:
                raise NameError('Server name not recognized.')

        # Call to system to submit the job array
        os.chdir(self.ParentFolder)     # Change into the directory so output files will be there
        os.system('qsub ' + os.path.join(self.ParentFolder, 'zacros_submit_JA.qs'))

        # Wait for jobs to be done
        all_jobs_done = False
        check_num = 1
        while not all_jobs_done:
            time.sleep(60)
            check_num += 1
            all_jobs_done = True
            for fldr in self.run_dirs:
                self.runtemplate.Path = fldr
                if not self.runtemplate.CheckComplete():
                    all_jobs_done = False

        sys.stdout.write('Jobs in ' + self.ParentFolder + ' have finished\n')
        sys.stdout.flush()


    def RunAllJobs_serial(self):       # Serial version of running all jobs

        '''
        Runs all trajectories in serial
        '''

        # loop over array of directories and execute the Zacros executable in that folder
        for fldr in self.run_dirs:
            self.runtemplate.Path = fldr
            self.runtemplate.Run_sim()


    def ReadMultipleRuns(self):

        '''
        Read all Zacros jobs in a given folder
        '''

        sys.stdout.write('Reading all runs in ' + self.ParentFolder + '\n')
        sys.stdout.flush()

        dummy_run = kmc_traj()       # Use this to transfer information


        if not self.run_dirs:   # If directory list is empty, fill it with the directories in the current folder which have finished jobs

            DirList = [d for d in os.listdir(self.ParentFolder) if os.path.isdir(os.path.join(self.ParentFolder, d))]      # List all folders in ParentFolder

            for direct in DirList:

                full_direct = os.path.join(self.ParentFolder, direct)
                dummy_run.Path = full_direct

                if dummy_run.CheckComplete():
                    self.run_dirs.append(full_direct)

            self.n_trajectories = len(self.run_dirs)

        # Create arrays for data
        # Use input data from runtemplate to properly size the arrays
        self.t_vec = []
        self.species_pops = []
        self.rxn_freqs = []
        self.History_final_snaps = []            # list of the final states
        self.propensities = []
        self.Props_integ = []
        self.traj_derivs = []
        self.events_total = []
        self.CPU_total = []
        self.time_avg_covs = []                   # surface species coverage based on time integral not including empty site
        #self.TS_site_props_sss = []

        for traj_dir in self.run_dirs:

            # Switch to folder and read output files
            dummy_run.Path = traj_dir
            dummy_run.ReadAllOutput()
            self.runtemplate = dummy_run            # so we have an example

            # Pass the data from this run into the larger data structures in
            # so you will not have to store separate kmc_traj objects
            # The large arrays of data will be easier to process

            self.t_vec = np.array(dummy_run.specnumout.t )
            self.species_pops.append(dummy_run.specnumout.spec )
            self.rxn_freqs.append(dummy_run.procstatout.events )

            if hasattr(dummy_run.histout, 'snapshots'):         # For some reason dummy_run.histout was not initializing properly...
                if not dummy_run.histout.snapshots == []:
                    self.History_final_snaps.append( dummy_run.histout.snapshots[-1] )

            self.propensities.append(dummy_run.prop)
            self.Props_integ.append(dummy_run.propCounter)
            self.traj_derivs.append(dummy_run.W_sen_anal)

            self.events_total.append( dummy_run.genout.events_occurred )
            self.CPU_total.append( dummy_run.genout.CPU_time )
            self.time_avg_covs.append(dummy_run.spec_num_int)
            #self.TS_site_props_sss.append(dummy_run.TS_site_props_ss)

        # Convert the data from lists to arrays
        self.species_pops = np.array(self.species_pops)
        self.rxn_freqs = np.array(self.rxn_freqs)
        self.propensities = np.array(self.propensities)
        self.Props_integ = np.array(self.Props_integ)
        self.traj_derivs = np.array(self.traj_derivs)
        self.events_total = np.array(self.events_total)
        self.CPU_total = np.array(self.CPU_total)
        self.time_avg_covs = np.array(self.time_avg_covs)
        #self.TS_site_props_sss = np.array(self.TS_site_props_sss)

        self.runAvg = copy.deepcopy(dummy_run)      # Initialize run average with information from dummyrun
        self.avg_updated = False


    def AverageRuns(self):

        '''
        Average multiple trajectories
        '''

        self.runAvg.Path = self.ParentFolder
        self.runAvg.specnumout.t = self.t_vec

        self.runAvg.specnumout.spec = np.mean(self.species_pops, axis = 0)
        self.runAvg.procstatout.events = np.mean(self.rxn_freqs, axis = 0)
        if not self.runAvg.prop is None:
            self.runAvg.prop = np.mean(self.propensities, axis = 0)
        if not self.runAvg.spec_num_int is None:
            self.runAvg.spec_num_int = np.mean(self.time_avg_covs,axis = 0)
        if not self.runAvg.propCounter is None:
            self.runAvg.propCounter = np.mean(self.Props_integ, axis = 0)
        #if not self.runAvg.TS_site_props_ss is None:
            #self.runAvg.TS_site_props_ss = np.mean(self.TS_site_props_sss, axis = 0)

        self.runAvg.genout.events_occurred = np.mean(self.events_total)
        self.runAvg.genout.CPU_time = np.mean(self.CPU_total)
        #self.runAvg.TS_site_props_ss = np.mean(self.TS_site_props_sss, axis = 0)

        self.avg_updated = True


    def Compute_batch_data(self, n_batches_total = 1000):

        '''
        Compute varaibles needed for the analysis, based on the input
        '''
        self.N_batches = n_batches_total
        # Find index of product molecule
        try:
            self.gas_prod_ind = len( self.runAvg.simin.surf_spec ) + self.runAvg.simin.gas_spec.index(self.gas_product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species ' + self.gas_product + ' not found.')

        # Find the stochiometry of the product molecule for each reaction
        nRxns = len(self.runAvg.genout.RxnNameList)
        self.TOF_stoich = np.zeros(nRxns)
        for i, elem_stoich in enumerate(self.runAvg.genout.Nu):
            self.TOF_stoich[i] = elem_stoich[self.gas_prod_ind]

        self.Nbpt = np.max( [ 3 , (self.N_batches-1) / self.n_trajectories + 2 ] )            # Set the number of batches per trajectory
        self.batch_length = self.t_vec[-1] / self.Nbpt


    def Compute_rate(self, include_ACF_CI = True):

        '''
        Use batch means to compute the reaction rate (and confidence interval)
        Also compute the autocorrelation function (ACF) (and confidence interval)

        :params include_ACF_CI: Whether to use statistical bootstrapping to compute confidence intervals.
            This takes some CPU time.
        :returns: The rate
        '''

        if not self.avg_updated:
            self.AverageRuns()

        self.Compute_batch_data()

        bin_edges = np.linspace(0, self.t_vec[-1], self.Nbpt + 1)
        rate_data = np.zeros([int(self.n_trajectories), int(self.Nbpt)])

        for i in range(self.Nbpt):

            idb_start = self.runAvg.time_search_interp(bin_edges[i])
            idb_end = self.runAvg.time_search_interp(bin_edges[i+1])

            prop_integ_start = idb_start[1][0] * self.Props_integ[:, idb_start[0][0], :] + idb_start[1][1] * self.Props_integ[:, idb_start[0][1], :]
            prop_integ_end = idb_end[1][0] * self.Props_integ[:, idb_end[0][0], :] + idb_end[1][1] * self.Props_integ[:, idb_end[0][1], :]

            rate_data[:,i] = np.dot ( ( prop_integ_end - prop_integ_start ) / self.batch_length , self.TOF_stoich )


        '''
        Compute confidence in the rate (assuming IID)
        '''

        all_batch_rates = rate_data[:,1::]
        all_batch_rates = all_batch_rates.flatten()

        confidence=0.90
        n = len(all_batch_rates)
        self.rate, se = np.mean(all_batch_rates), scipy.stats.sem(all_batch_rates)
        self.rate_CI = se * scipy.stats.t._ppf((1+confidence)/2., n - 1)

        '''
        Compute autocorrelation function
        '''

        # Compute the autocorrelation function for the batch means of the rate
        c1 = rate_data[:,1:-1:]
        c1 = c1.flatten()
        c2 = rate_data[:,2::]
        c2 = c2.flatten()


        if np.var(all_batch_rates) == 0:        # If the variance of the batch means is zero, autocorrelation cannot be computed.
            self.ACF = None
            self.ACF_CI = None
        else:

            self.ACF = ( np.mean(c1 * c2) - np.mean(c1) * np.mean(c2) ) / np.var(all_batch_rates)

            if include_ACF_CI:

                '''
                Compute confidence interval for ACF
                '''

                N_boot = 100

                ACF_dist = np.zeros(N_boot)
                for i in range(N_boot):
                    subpop_inds = np.random.randint(len(c1), size=len(c1))
                    c1_new = c1[subpop_inds]
                    c2_new = c2[subpop_inds]
                    ACF_dist[i] = ( np.mean(c1_new * c2_new) - np.mean(c1_new) * np.mean(c2_new) ) / np.var(all_batch_rates)

                ACF_dist = sorted(ACF_dist)
                ACF_high = ACF_dist[int(0.95 * N_boot)]
                ACF_low = ACF_dist[int(0.05 * N_boot)]
                self.ACF_CI = (ACF_high - ACF_low) / 2

        # Rescale error bars on the rate based on the ACF
        if not self.ACF is None:
            var_scale = 1
            for i in range(1,self.Nbpt):
                var_scale = var_scale + 2 * ( 1 - i / ( self.Nbpt - 1 ) ) * self.ACF ** i
            self.rate_CI = self.rate_CI * np.sqrt( var_scale )

        return self.rate


    def PerformSA(self, delta_t = None, ergodic = True, dp_per_bin = 10):          # Need implement time point interpolation

        '''
        Perform likelihood ratio sensitivity analysis with a combination of time and trajectory averaging
        :param delta_t:   Size of the time window used for likelihood ratio sensitivity analysis. By default, it is the size of a batch.
        :param ergodic:   True - average the rate over the entire time interval (centered ergodic likelihood ratio)
                    False - use the rate at the end of the time interval (centered likelihood ratio)
         Data between sample points is estimated with linear interpolation
        '''

        self.Compute_batch_data()

        # Use entire trajectory length as the default time window
        if delta_t is None:
            delta_t = self.batch_length

        if not self.avg_updated:
            self.AverageRuns()

        if delta_t > self.t_vec[-1] or delta_t < 0:
            print('Time window: ' + str(delta_t))
            print('Final time: ' + str(self.t_vec[-1]))
            raise Exception('Time window is too large. Insufficient sampling.')


        bin_edges = np.linspace(0, self.t_vec[-1], self.Nbpt * dp_per_bin + 1)

        dp_per_traj = dp_per_bin * (self.Nbpt-1) + 1
        rate_data_erg = np.zeros( self.n_trajectories * dp_per_traj  )
        W_data_all = np.zeros([ self.traj_derivs.shape[2] , self.n_trajectories * dp_per_traj ] )
        fW_data = np.zeros( self.n_trajectories * dp_per_traj )
        rate_contributions_all = np.zeros(self.traj_derivs.shape[2])


        data_ind = 0
        for traj_ind in range(self.n_trajectories):
            for i in range(dp_per_traj):

                idb_start = self.runAvg.time_search_interp( bin_edges[i] )
                idb_end = self.runAvg.time_search_interp( bin_edges[i+dp_per_bin] )

                W_start = idb_start[1][0] * self.traj_derivs[ traj_ind, idb_start[0][0] , :] + idb_start[1][1] * self.traj_derivs[ traj_ind, idb_start[0][1] , :]
                W_end = idb_end[1][0] * self.traj_derivs[ traj_ind, idb_end[0][0] , :] + idb_end[1][1] * self.traj_derivs[ traj_ind, idb_end[0][1] , :]
                W_data_all[:, data_ind] = W_end - W_start

                if not ergodic:
                    idb_start = self.runAvg.time_search_interp( bin_edges[i+dp_per_bin-1] )

                if ergodic:
                    prop_integ_start = idb_start[1][0] * self.Props_integ[traj_ind, idb_start[0][0], :] + idb_start[1][1] * self.Props_integ[traj_ind, idb_start[0][1], :]
                    prop_integ_end = idb_end[1][0] * self.Props_integ[traj_ind, idb_end[0][0], :] + idb_end[1][1] * self.Props_integ[traj_ind, idb_end[0][1], :]
                    rate_data_erg[data_ind] = np.dot ( ( prop_integ_end - prop_integ_start ) / self.batch_length , self.TOF_stoich )
                    rate_contributions_all = rate_contributions_all + ( ( prop_integ_end - prop_integ_start ) / self.batch_length * self.TOF_stoich )
                else:
                    inst_rates = idb_end[1][0] * self.propensities[traj_ind, idb_end[0][0], :] + idb_end[1][1] * self.propensities[traj_ind, idb_end[0][1], :]
                    rate_data_erg[data_ind] = np.dot ( inst_rates , self.TOF_stoich )
                    rate_contributions_all = rate_contributions_all + ( inst_rates * self.TOF_stoich )

                data_ind += 1

        # Normalize rate contributions by the number of data points
        rate_contributions_all = rate_contributions_all / ( self.n_trajectories * dp_per_traj )

        # Combine forward and reverse reactions
        W_data = np.zeros([self.runAvg.mechin.get_num_rxns(), self.n_trajectories * dp_per_traj])
        rate_contributions = np.zeros(self.runAvg.mechin.get_num_rxns())
        ind = 0
        for i in range(self.runAvg.mechin.get_num_rxns()):

            rxn_and_var = self.runAvg.mechin.get_rxn_var_inds(i)

            if self.runAvg.mechin.rxn_list[rxn_and_var[0]].is_reversible:
                W_data[i, :] = W_data_all[ind, :] + W_data_all[ind+1, :]
                rate_contributions[i] = rate_contributions_all[ind] + rate_contributions_all[ind+1]
                ind += 2
            else:
                W_data[i, :] = W_data_all[ind, :]
                rate_contributions[i] = rate_contributions_all[i]
                ind += 1

        # Calculate NSCs
        mean_rate = 0
        W_mean = np.zeros( self.runAvg.mechin.get_num_rxns() )
        NSCs = np.zeros(self.runAvg.mechin.get_num_rxns())

        for dp in range( self.n_trajectories * dp_per_traj ):
            mean_rate = mean_rate + rate_data_erg[dp]
            W_mean = W_mean + W_data[ : , dp ]
            NSCs = NSCs + rate_data_erg[dp] * W_data[ : , dp ]

        # Normalize means
        mean_rate = mean_rate / ( self.n_trajectories * dp_per_traj )
        W_mean = W_mean / ( self.n_trajectories * dp_per_traj )
        NSCs = NSCs / ( self.n_trajectories * dp_per_traj )

        NSCs = NSCs - W_mean * mean_rate + rate_contributions       # Convert from ELR to CELR
        NSCs = NSCs / mean_rate     # normalize


        '''
        Compute error bounds on NSCs
        '''

        N_boot = 100

        NSC_sam = np.zeros([ N_boot , len(NSCs) ])
        for i in range(N_boot):

            # Prepare new set of data
            subpop_inds = np.random.randint(self.n_trajectories * dp_per_traj, size = self.n_trajectories * dp_per_traj )

            rate_data_erg_sub = rate_data_erg[subpop_inds]
            W_sub = W_data[:, subpop_inds]

            # Calculate NSCs
            mean_rate = 0
            W_mean = np.zeros(self.runAvg.mechin.get_num_rxns())
            NSCsub = np.zeros(self.runAvg.mechin.get_num_rxns())

            for dp in range( self.n_trajectories * dp_per_traj ):
                mean_rate = mean_rate + rate_data_erg_sub[dp]
                W_mean = W_mean + W_sub[ : , dp ]
                NSCsub = NSCsub + rate_data_erg_sub[dp] * W_sub[ : , dp ]

            # Normalize means
            mean_rate = mean_rate / ( self.n_trajectories * dp_per_traj )
            W_mean = W_mean / ( self.n_trajectories * dp_per_traj )
            NSCsub = NSCsub / ( self.n_trajectories * dp_per_traj )

            NSCsub = NSCsub - W_mean * mean_rate + rate_contributions       # Convert from ELR to CELR
            NSCsub = NSCsub / mean_rate     # normalize

            NSC_sam[i,:] = NSCsub

        # Sort each column of the data and compute confidence interval
        NSC_ci = np.zeros(NSCs.shape)
        for rxn_ind in range(len(NSCs)):
            NSC_dist = sorted(NSC_sam[:,rxn_ind])
            NSC_dist_high = NSC_dist[int(0.95 * N_boot)]
            NSC_dist_low = NSC_dist[int(0.05 * N_boot)]
            NSC_ci[rxn_ind] = (NSC_dist_high - NSC_dist_low) / 2

        self.NSC = NSCs
        self.NSC_ci = NSC_ci


    def PlotSensitivities(self, NSC_cut = 0.05):

        '''
        Plot the results of sensitivity analysis
        NSC_cut : Reactions with normalized sensitivity coefficients (NSC) below this limit will not be plotted
        '''

        PlotOptions()
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []

        for i in range (self.runAvg.mechin.get_num_rxns()):

            if self.NSC[i] + self.NSC_ci[i] > NSC_cut or self.NSC[i] - self.NSC_ci[i] < -NSC_cut:
                plt.barh(ind-0.9, self.NSC[i], width, color='r', xerr = self.NSC_ci[i], ecolor='k')

                indices = get_rxn_var_inds(i)
                ylabels.append(self.runAvg.rxn_list[indices[0]].name + '_' + self.runAvg.rxn_list[indices[1]].name )
                yvals.append(ind-0.6)
                ind = ind - 1

        plt.plot([0, 0], [0, ind], color='k')
        plt.xlim([0,1])
        plt.xticks(size=20)
        plt.yticks(size=10)
        plt.xlabel('NSC',size=24)
        plt.yticks(yvals, ylabels)
        ax = plt.subplot(111)
        pos = [0.30, 0.15, 0.65, 0.8]
        ax.set_position(pos)

        plt.savefig(os.path.join(self.ParentFolder, 'SA_output.png'))
        plt.close()


    def WriteSA_output(self):

        '''
        Write the results of sensitivity analysis into a SA_output.txt
        '''

        with open(os.path.join(self.ParentFolder, 'SA_output.txt'), 'w') as txt:
            txt.write('Normalized sensitivity coefficients \n\n')
            #txt.write('Turnover frequency: \t' + '{0:.3E} \t'.format(self.TOF) + '+- {0:.3E} \t'.format(self.TOF_error) + '\n\n')
            txt.write('Reaction name \t NSC \t NSC confidence \n')

            ind = 0
            for rxn in self.runAvg.mechin.rxn_list:
                for vrnt in rxn.variant_list:
                    txt.write(rxn.name + '_' + vrnt.name + '\t' + '{0:.3f} +- \t'.format(self.NSC[ind]) + '\n')


def append_replicates(batch1, batch2):

    '''
    Append two sets of trajectories
    :returns: Replicates object with all trajectories appended
    '''

    sand = copy.deepcopy(batch2)

    sand.t_vec = np.hstack( [ batch1.t_vec, batch2.t_vec[1::] + batch1.t_vec[-1]] )

    #sand.History_final_snaps = batch2.History_final_snaps          # history will already have been copied
    sand.events_total = batch1.events_total + batch2.events_total
    sand.CPU_total = batch1.CPU_total + batch2.CPU_total

    # Cumulative data structures

    # Add data from the end of the first calculation to the second calculation
    for traj_ind in range(batch2.n_trajectories):
        for time_ind in range(len(batch2.t_vec)):

            for rxn_ind in range( sand.rxn_freqs.shape[2] ):

                sand.rxn_freqs[traj_ind, time_ind, rxn_ind] += batch1.rxn_freqs[traj_ind, -1, rxn_ind]
                sand.Props_integ[traj_ind, time_ind, rxn_ind] += batch1.Props_integ[traj_ind, -1, rxn_ind]
                sand.traj_derivs[traj_ind, time_ind, rxn_ind] += batch1.traj_derivs[traj_ind, -1, rxn_ind]

            # Add to the end for gas phase populations
            for spec_ind in range( len( batch2.runAvg.simin.surf_spec ) , len( batch2.runAvg.simin.surf_spec ) + len( batch2.runAvg.simin.gas_spec ) ):
                sand.species_pops[traj_ind, time_ind, spec_ind] += batch1.species_pops[traj_ind, -1, spec_ind]

    # Combine the data
    sand.species_pops = np.concatenate( [batch1.species_pops, sand.species_pops[:,1::,:]], axis = 1 )
    sand.rxn_freqs = np.concatenate( [batch1.rxn_freqs, sand.rxn_freqs[:,1::,:]], axis = 1 )
    sand.propensities = np.concatenate( [batch1.propensities, sand.propensities[:,1::,:]], axis = 1 )
    sand.Props_integ = np.concatenate( [batch1.Props_integ, sand.Props_integ[:,1::,:]], axis = 1 )
    sand.traj_derivs = np.concatenate( [batch1.traj_derivs, sand.traj_derivs[:,1::,:]], axis = 1 )
    #sand.TS_site_props_sss = ( batch1.TS_site_props_sss * batch1.t_vec[-1] + batch2.TS_site_props_sss * batch2.t_vec[-1] ) / (batch1.t_vec[-1] + batch2.t_vec[-1])

    sand.avg_updated = False

    return sand
