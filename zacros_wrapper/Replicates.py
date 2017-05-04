import os
import numpy as np
import matplotlib as mat
#mat.use('Agg')
import matplotlib.pyplot as plt
import copy
import sys
import scipy.stats      # causes trouble with newer version of Python

from KMC_Run import kmc_traj
from Helper import *
import time
import scipy


class Replicates:

    '''
    Has data from muliple trajectories with the same input, but different 
    random seeds and possibly different initial states
    This class handles all of the statistics
    '''
    
    def __init__(self):
             
        # General info
        self.ParentFolder = None
        self.runtemplate = None                             # Use this to build replicate jobs
        self.runAvg = None            # Values are averages of all runs from runList
        self.avg_updated = False        # Whether we need to update the trajectory averages
        self.n_runs = 0
        self.N_batches = 1000       # Used to be 1000
        self.Nbpt = None
        self.batch_length = None    # Length of time in one batch
        self.gas_product = None
        
        # Input data for different trajectories
        self.run_dirs = []
        self.rand_seeds = []
        
        # Output data taken accross all trajectories
        # Put here various multidimensional numpy arrays which will make it easier to do statistics
        self.t_vec = None                    # times at which data is recorded
        self.species_pops = None
        self.rxn_freqs = None
        self.History_final_snaps = None
        #self.props = None            # We do not print out this output file
        self.Props_integ = None
        self.traj_derivs = None
        self.events_total = None
        self.CPU_total = None
        
        # Needed for computing the rates
        self.gas_prod_ind = None
        self.TOF_stoich = None
        
        # Analysis
        self.rates_data = None
        self.TOF = None
        self.TOF_error = None
        self.NSC = None
        self.NSC_ci = None
    
    
    def Set_kmc_template(self, kmc_temp):
        
        '''
        Set the KMC trajectory template as well as other useful variables
        '''
        
        self.runtemplate = copy.deepcopy(kmc_temp)
        self.gas_product = kmc_temp.gas_prod
        
    
    def BuildJobFiles(self, init_states = []):
        
        '''
        Builds folders with Zacros input files. The only difference is the random seed.
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
        for i in range(self.n_runs):
            
            # Set the random seed
            self.runtemplate.Conditions['Seed'] = seed
            self.rand_seeds.append(seed)
            seed = seed + 1
            
            # Set the path
            fldr = os.path.join(self.ParentFolder, str(i+1))
            self.runtemplate.Path = fldr
            self.run_dirs.append(fldr)
            
            # Set initial state
            if not init_states == []:
                self.runtemplate.StateInput['Type'] = 'history'
                self.runtemplate.StateInput['Struct'] = init_states[i]
                
            # Write input files in that folder
            if not os.path.exists(fldr):
                os.makedirs(fldr)
                
            self.runtemplate.WriteAllInput()
            
        
    def RunAllJobs_parallel_JobArray(self, max_cores = 100, server = 'Squidward'):
    
        '''
        Runs a job array on Squidward or Farber
        Writes the .qs file with certain parameters changed
        Then, it submits it to the gridengine and waits for them to finish
        '''
    
        sys.stdout.write('Running parallel jobs')
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
        
        n_cores = np.min([max_cores, self.n_runs])           
        
        with open(os.path.join(self.ParentFolder, 'zacros_submit_JA.qs'), 'w') as txt:
            
            if server == 'Farber':             
            
                txt.write('#!/bin/bash\n')
                txt.write('#$ -cwd\n')
                txt.write('#$ -j y\n')
                txt.write('#$ -S /bin/bash\n')
                txt.write('#$ -l h_cpu=168:00:00\n')
                txt.write('#\n')
                txt.write('\n')
                txt.write('#$ -N zacros_JA 					#This is the name of the job array\n')
                txt.write('#$ -t 1-' + str(self.n_runs) + '  							#Assumes task IDs increment by 1; can also increment by another value\n')
                txt.write('#$ -tc ' + str(n_cores) + ' 							#This is the total number of tasks to run at any given moment\n')
                txt.write('#$ -pe threads 1 				#Change the last field to the number of processors desired per task\n')
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
            
            else:       # Squidward
            
                txt.write('#!/bin/bash\n')
                txt.write('#$ -cwd\n')
                txt.write('#$ -j y\n')
                txt.write('#$ -S /bin/bash\n')
                txt.write('#\n')
                txt.write('\n')
                txt.write('#$ -N zacros_JA 					#This is the name of the job array\n')
                txt.write('#$ -t 1-' + str(self.n_runs) + '  							#Assumes task IDs increment by 1; can also increment by another value\n')
                txt.write('#$ -tc ' + str(n_cores) + ' 							#This is the total number of tasks to run at any given moment\n')
                txt.write('#$ -pe openmpi-smp 1 				#Change the last field to the number of processors desired per task\n')
                txt.write('\n')
                txt.write('job_file=\'' + os.path.join(self.ParentFolder, 'dir_list.txt') + '\'\n')
                txt.write('#Change to the job directory\n')
                txt.write('job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")\n')
                txt.write('cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>\n')
                txt.write('\n\n')
                txt.write('time ' + self.runtemplate.exe_file)
    
        # Call to system to submit the job array
        os.chdir(self.ParentFolder)
        os.system('qsub ' + os.path.join(self.ParentFolder, 'zacros_submit_JA.qs'))
    
        # Wait for jobs to be done
        all_jobs_done = False
        while not all_jobs_done:
            time.sleep(60)
            all_jobs_done = True
            for fldr in self.run_dirs:
                self.runtemplate.Path = fldr
                if not self.runtemplate.CheckComplete():
                    all_jobs_done = False
                    
        sys.stdout.write('Jobs in ' + self.ParentFolder + ' have finished\n')
        sys.stdout.flush()
        
    
    def RunAllJobs_serial(self):       # Serial version of running all jobs

        '''
        Runs all Zacros jobs in serial
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
                    
            self.n_runs = len(self.run_dirs)
        
        # Create arrays for data
        # Use input data from runtemplate to properly size the arrays
        self.t_vec = []
        self.species_pops = []
        self.rxn_freqs = []
        self.History_final_snaps = []            # list of the final states
        #self.props = []
        self.Props_integ = []
        self.traj_derivs = []
        self.events_total = []
        self.CPU_total = []
        
        for i in range(self.n_runs):
        
            # Switch to folder and read output files
            dummy_run.Path = self.run_dirs[i]
            dummy_run.ReadAllOutput()
            self.runtemplate = dummy_run            # so we have an example
            
            # Pass the data from this run into the larger data structures in 
            # so you will not have to store separate kmc_traj objects
            # The large arrays of data will be easier to process
            
            self.t_vec = np.array(dummy_run.Specnum['t'])
            self.species_pops.append(dummy_run.Specnum['spec'])
            self.rxn_freqs.append(dummy_run.Procstat['events'])
            
            if not dummy_run.History == []:
                self.History_final_snaps.append(dummy_run.History[-1])
            
            #self.props.append(dummy_run.Binary['prop'])
            self.Props_integ.append(dummy_run.Binary['propCounter'])
            self.traj_derivs.append(dummy_run.Binary['W_sen_anal'])
            
            self.events_total.append(dummy_run.Performance['events_occurred'])
            self.CPU_total.append(dummy_run.Performance['CPU_time'])
        
        # Convert the data from lists to arrays
        self.species_pops = np.array(self.species_pops)
        self.rxn_freqs = np.array(self.rxn_freqs)
        #self.props = np.array(self.props)
        self.Props_integ = np.array(self.Props_integ)
        self.traj_derivs = np.array(self.traj_derivs)
        self.events_total = np.array(self.events_total)
        self.CPU_total = np.array(self.CPU_total)
        
        self.runAvg = copy.deepcopy(dummy_run)      # Initialize run average with information from dummyrun
        self.avg_updated = False
        

    def Set_analysis_varaibles(self):
    
        '''
        Compute varaibles needed for the analysis, based on the input
        '''
    
        # Find index of product molecule
        try:
            self.gas_prod_ind = self.runAvg.Species['n_surf'] + self.runAvg.Species['gas_spec'].index(self.gas_product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species ' + self.gas_product + ' not found.')
        
        # Find the stochiometry of the product molecule for each reaction         
        nRxns = len(self.runAvg.Reactions['Nu'])
        self.TOF_stoich = np.zeros(nRxns)
        for i, elem_stoich in enumerate(self.runAvg.Reactions['Nu']):
            self.TOF_stoich[i] = elem_stoich[self.gas_prod_ind]
            
        self.Nbpt = np.max( [ 3 , (self.N_batches-1) / self.n_runs + 2 ] )            # Set the number of batches per trajectory
        self.batch_length = self.t_vec[-1] / self.Nbpt
        
        
    def AverageRuns(self):
        
        '''
        Average multiple trajectories
        '''
        
        self.runAvg.Path = self.ParentFolder
        self.runAvg.Specnum['t'] = self.t_vec
        
        self.runAvg.Specnum['spec'] = np.mean(self.species_pops, axis = 0)
        self.runAvg.Procstat['events'] = np.mean(self.rxn_freqs, axis = 0)
        
        #self.runAvg.Binary['prop'] = np.mean(self.props, axis = 0)
        if not self.runAvg.Binary['propCounter'] is None:
            self.runAvg.Binary['propCounter'] = np.mean(self.Props_integ, axis = 0)
        
        self.runAvg.Performance['events_occurred'] = np.mean(self.events_total)
        self.runAvg.Performance['CPU_time'] = np.mean(self.CPU_total)
        
        self.avg_updated = True
        
    
    def Compute_ACF(self):
    
        '''
        Check steady state on the basis of batch means
        '''
        
        if not self.avg_updated:
            self.AverageRuns()
        
        self.Set_analysis_varaibles()
            
        bin_edges = np.linspace(0, self.t_vec[-1], self.Nbpt + 1)
        rate_data = np.zeros([self.n_runs, self.Nbpt])
        
        for i in range(self.Nbpt):
        
            idb_start = self.runAvg.time_search_interp(bin_edges[i])
            idb_end = self.runAvg.time_search_interp(bin_edges[i+1])
            
            prop_integ_start = idb_start[1][0] * self.Props_integ[:, idb_start[0][0], :] + idb_start[1][1] * self.Props_integ[:, idb_start[0][1], :]
            prop_integ_end = idb_end[1][0] * self.Props_integ[:, idb_end[0][0], :] + idb_end[1][1] * self.Props_integ[:, idb_end[0][1], :]
            
            rate_data[:,i] = np.dot ( ( prop_integ_end - prop_integ_start ) / self.batch_length , self.TOF_stoich )
        
        '''
        Plot Rates
        '''
        
        PlotOptions()
        plt.figure()
        
        for i in range(self.n_runs):
            plt.plot(bin_edges[2::], rate_data[i,1::], '.')
        
        #plt.xticks(size=20)
        #plt.yticks(size=20)
        plt.xlabel('Time (s)', size=24)
        plt.ylabel('Rate', size=24)
        
        ax = plt.subplot(111)
        ax.set_position([0.2, 0.15, 0.7, 0.8])
        
        
        plt.savefig(os.path.join(self.ParentFolder, 'bin_rates.png'))
        plt.close()
        
        
        plt.figure()
        
        plt.plot(bin_edges[2::], np.mean( rate_data[:,1::], axis=0 ) , '.')
        
        #plt.xticks(size=20)
        #plt.yticks(size=20)
        plt.xlabel('Time (s)', size=24)
        plt.ylabel('Rate', size=24)
        
        ax = plt.subplot(111)
        ax.set_position([0.2, 0.15, 0.7, 0.8])
        
        
        plt.savefig(os.path.join(self.ParentFolder, 'mean_bin_rates.png'))
        plt.close()
        
        '''
        Compute confidence in the rate (assuming IID)
        '''
        
        alldata = rate_data[:,1::]
        alldata = alldata.reshape(np.prod(alldata.shape))

        confidence=0.90
        n = len(alldata)
        mean_rate, se = np.mean(alldata), scipy.stats.sem(alldata)
        rate_ci = se * scipy.stats.t._ppf((1+confidence)/2., n - 1)
        
        '''
        Compute ACF
        '''
        
        # Compute the autocorrelation function for the batch means of the rate
        c1 = rate_data[:,1:-1:]
        c1 = c1.reshape( np.prod(c1.shape) )
        c2 = rate_data[:,2::]
        c2 = c2.reshape( np.prod(c2.shape) )
        
        
        if ( np.var(alldata) == 0 ) or ( np.mean(rate_data[ : , 1 ]) == 0 ):        # Covers the case where we have zero rate
            ACF = None
            ACF_ci = None
        else:
            ACF = ( np.mean(c1 * c2) - np.mean(c1) * np.mean(c2) ) / np.var(alldata)


            '''
            Compute error bounds on ACF
            '''
            
            N_boot = 100
            
            ACF_dist = np.zeros(N_boot)
            for i in range(N_boot):
                subpop_inds = np.random.randint(len(c1), size=len(c1))
                c1_new = c1[subpop_inds]
                c2_new = c2[subpop_inds]
                ACF_dist[i] = ( np.mean(c1_new * c2_new) - np.mean(c1_new) * np.mean(c2_new) ) / np.var(alldata)
                #ACF_dist[i] = ( np.mean(c1_new * c2_new) - np.mean(alldata) ** 2 ) / np.var(alldata)
                #ACF_dist[i] = np.mean(c1_new * c2_new) - np.mean(c1_new) * np.mean(c2_new)              # test not normalizing the rate
                
            ACF_dist = sorted(ACF_dist)
            ACF_high = ACF_dist[int(0.95 * N_boot)]
            ACF_low = ACF_dist[int(0.05 * N_boot)]
            ACF_ci = (ACF_high - ACF_low) / 2
        
        return {'Rate': mean_rate, 'Rate_ci': rate_ci, 'ACF': ACF, 'ACF_ci': ACF_ci}
        #return [ACF, ACF_high, ACF_low]
        

    def Compute_rate(self, delta_t = None):          # Need implement time point interpolation
        
        '''
        Perform sensitivity analysis with a combination of
        time and trajectory averaging
        
        delta_t             size of the time window
        ergodic             True: use the rate at the end of the time interval, False: average the rate over the entire time interval 
        
        The time-averaging could be improved in here. Right now it is a little weird...
        You can use linear interpolation to get data between time points
        Can also make it so that the statistics make sense for non-evenly spaced time points
        '''
        
        self.Set_analysis_varaibles()
        
        # Use entire trajectory length as the default time window
        if delta_t is None:
            delta_t = self.batch_length
        
        if not self.avg_updated:
            self.AverageRuns()
            
        if delta_t > self.t_vec[-1] or delta_t < 0:
            print 'Time window: ' + str(delta_t)
            print 'Final time: ' + str(self.t_vec[-1])
            raise Exception('Time window is too large. Insufficient sampling.')
        
        bin_edges = np.linspace(0, self.t_vec[-1], self.Nbpt + 1)
        
        rate_data_erg = np.zeros( self.n_runs )
        rate_contributions_all = np.zeros(self.traj_derivs.shape[2])
         
        
        data_ind = 0
        for traj_ind in range(self.n_runs):
            
            idb_start = self.runAvg.time_search_interp( bin_edges[1] )
            idb_end = self.runAvg.time_search_interp( bin_edges[-1] )
            
            prop_integ_start = idb_start[1][0] * self.Props_integ[traj_ind, idb_start[0][0], :] + idb_start[1][1] * self.Props_integ[traj_ind, idb_start[0][1], :]
            prop_integ_end = idb_end[1][0] * self.Props_integ[traj_ind, idb_end[0][0], :] + idb_end[1][1] * self.Props_integ[traj_ind, idb_end[0][1], :]
            rate_data_erg[data_ind] = np.dot ( ( prop_integ_end - prop_integ_start ) / ( bin_edges[-1] - bin_edges[1] ) , self.TOF_stoich )
          
            data_ind += 1
        
        # Calculate NSCs
        mean_rate = 0
        
        for dp in range( self.n_runs):
            mean_rate = mean_rate + rate_data_erg[dp]
            
        # Normalize means
        mean_rate = mean_rate / self.n_runs
        
        return mean_rate

        
    def PerformSA(self, delta_t = None, ergodic = True):          # Need implement time point interpolation
        
        '''
        Perform sensitivity analysis with a combination of
        time and trajectory averaging
        
        delta_t             size of the time window
        ergodic             True: use the rate at the end of the time interval, False: average the rate over the entire time interval 
        
        The time-averaging could be improved in here. Right now it is a little weird...
        You can use linear interpolation to get data between time points
        Can also make it so that the statistics make sense for non-evenly spaced time points
        '''
        
        self.Set_analysis_varaibles()
        
        # Use entire trajectory length as the default time window
        if delta_t is None:
            delta_t = self.batch_length
        
        if not self.avg_updated:
            self.AverageRuns()
            
        if delta_t > self.t_vec[-1] or delta_t < 0:
            print 'Time window: ' + str(delta_t)
            print 'Final time: ' + str(self.t_vec[-1])
            raise Exception('Time window is too large. Insufficient sampling.')
        
        dp_per_bin = 10
        bin_edges = np.linspace(0, self.t_vec[-1], self.Nbpt * dp_per_bin + 1)
        
        dp_per_traj = dp_per_bin * (self.Nbpt-1) + 1 
        rate_data_erg = np.zeros( self.n_runs * dp_per_traj  )
        W_data_all = np.zeros([ self.traj_derivs.shape[2] , self.n_runs * dp_per_traj ] )
        fW_data = np.zeros( self.n_runs * dp_per_traj )
        rate_contributions_all = np.zeros(self.traj_derivs.shape[2])
         
        
        data_ind = 0
        for traj_ind in range(self.n_runs):
            for i in range(dp_per_traj):
            
                idb_start = self.runAvg.time_search_interp( bin_edges[i] )
                idb_end = self.runAvg.time_search_interp( bin_edges[i+dp_per_bin] )
                
                W_start = idb_start[1][0] * self.traj_derivs[ traj_ind, idb_start[0][0] , :] + idb_start[1][1] * self.traj_derivs[ traj_ind, idb_start[0][1] , :]
                W_end = idb_end[1][0] * self.traj_derivs[ traj_ind, idb_end[0][0] , :] + idb_end[1][1] * self.traj_derivs[ traj_ind, idb_end[0][1] , :]
                W_data_all[:, data_ind] = W_end - W_start
                
                if not ergodic:
                    idb_start = self.runAvg.time_search_interp( bin_edges[i+dp_per_bin-1] )
                
                prop_integ_start = idb_start[1][0] * self.Props_integ[traj_ind, idb_start[0][0], :] + idb_start[1][1] * self.Props_integ[traj_ind, idb_start[0][1], :]
                prop_integ_end = idb_end[1][0] * self.Props_integ[traj_ind, idb_end[0][0], :] + idb_end[1][1] * self.Props_integ[traj_ind, idb_end[0][1], :]
                rate_data_erg[data_ind] = np.dot ( ( prop_integ_end - prop_integ_start ) / self.batch_length , self.TOF_stoich )
                rate_contributions_all = rate_contributions_all + ( ( prop_integ_end - prop_integ_start ) / self.batch_length * self.TOF_stoich )
            
                data_ind += 1
        
        # Normalize rate contributions by the number of data points
        rate_contributions_all = rate_contributions_all / ( self.n_runs * dp_per_traj )
        
        # Combine forward and reverse reactions
        W_data = np.zeros([self.runAvg.Reactions['nrxns'], self.n_runs * dp_per_traj])
        rate_contributions = np.zeros(self.runAvg.Reactions['nrxns'])
        ind = 0
        for i in range(self.runAvg.Reactions['nrxns']):
            
            if self.runAvg.Reactions['is_reversible'][i]:
                W_data[i, :] = W_data_all[ind, :] + W_data_all[ind+1, :]
                rate_contributions[i] = rate_contributions_all[ind] + rate_contributions_all[ind+1]
                ind += 2
            else:
                W_data[i, :] = W_data_all[ind, :]
                rate_contributions[i] = rate_contributions_all[i]
                ind += 1
        
        # Calculate NSCs
        mean_rate = 0
        W_mean = np.zeros(self.runAvg.Reactions['nrxns'])
        NSCs = np.zeros(self.runAvg.Reactions['nrxns'])
        
        for dp in range( self.n_runs * dp_per_traj ):
            mean_rate = mean_rate + rate_data_erg[dp]
            W_mean = W_mean + W_data[ : , dp ]
            NSCs = NSCs + rate_data_erg[dp] * W_data[ : , dp ]
            
        # Normalize means
        mean_rate = mean_rate / ( self.n_runs * dp_per_traj )
        W_mean = W_mean / ( self.n_runs * dp_per_traj )
        NSCs = NSCs / ( self.n_runs * dp_per_traj )
        
        NSCs = NSCs - W_mean * mean_rate + rate_contributions       # Convert from ELR to CELR
        NSCs = NSCs / mean_rate     # normalize 
        
        print NSCs
        
        '''
        Compute error bounds on NSCs
        '''
        
        N_boot = 100
        
        NSC_sam = np.zeros([ N_boot , len(NSCs) ])
        for i in range(N_boot):
        
            # Prepare new set of data
            subpop_inds = np.random.randint(self.n_runs * dp_per_traj, size = self.n_runs * dp_per_traj )
            
            rate_data_erg_sub = rate_data_erg[subpop_inds]
            W_sub = W_data[:, subpop_inds]
            
            # Calculate NSCs
            mean_rate = 0
            W_mean = np.zeros(self.runAvg.Reactions['nrxns'])
            NSCsub = np.zeros(self.runAvg.Reactions['nrxns'])
            
            for dp in range( self.n_runs * dp_per_traj ):
                mean_rate = mean_rate + rate_data_erg_sub[dp]
                W_mean = W_mean + W_sub[ : , dp ]
                NSCsub = NSCsub + rate_data_erg_sub[dp] * W_sub[ : , dp ]
                
            # Normalize means
            mean_rate = mean_rate / ( self.n_runs * dp_per_traj )
            W_mean = W_mean / ( self.n_runs * dp_per_traj )
            NSCsub = NSCsub / ( self.n_runs * dp_per_traj )
            
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
        
        print NSC_ci
        
        
    def PlotSensitivities(self, NSC_cut = 0.05): 
        
        '''
        Plot the results of sensitivity analysis
        '''
        
        PlotOptions()
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        
        for i in range (self.runAvg.Reactions['nrxns']):
            
            if self.NSC[i] + self.NSC_ci[i] > NSC_cut or self.NSC[i] - self.NSC_ci[i] < -NSC_cut:     
                plt.barh(ind-0.9, self.NSC[i], width, color='r', xerr = self.NSC_ci[i], ecolor='k')               
                ylabels.append(self.runAvg.Reactions['names'][i])              
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
        Write the results of sensitivity analysis into an SA_output.txt
        '''
    
        with open(os.path.join(self.ParentFolder, 'SA_output.txt'), 'w') as txt:
            txt.write('Normalized sensitivity coefficients \n\n')
            #txt.write('Turnover frequency: \t' + '{0:.3E} \t'.format(self.TOF) + '+- {0:.3E} \t'.format(self.TOF_error) + '\n\n')               
            txt.write('Reaction name \t NSC \t NSC confidence \n')

            for rxn_ind in range(self.runAvg.Reactions['nrxns']):
                txt.write(self.runAvg.Reactions['names'][rxn_ind] + '\t' + '{0:.3f} +- \t'.format(self.NSC[rxn_ind]) + '\n')
            
            
    @staticmethod
    def time_sandwich(batch1, batch2):
        
        '''
        Append two sets of trajectories
        '''
        
        sand = copy.deepcopy(batch2)
        
        sand.t_vec = np.hstack( [ batch1.t_vec, batch2.t_vec[1::] + batch1.t_vec[-1]] ) 
        
        #sand.History_final_snaps = batch2.History_final_snaps          # history will already have been copied
        sand.events_total = batch1.events_total + batch2.events_total
        sand.CPU_total = batch1.CPU_total + batch2.CPU_total
        
        # Cumulative data structures
        
        # Add data from the end of the first calculation to the second calculation
        for traj_ind in range(batch2.n_runs):
            for time_ind in range(len(batch2.t_vec)):
            
                for rxn_ind in range(len(batch2.runAvg.Reactions['Nu'])):
            
                    sand.rxn_freqs[traj_ind, time_ind, rxn_ind] += batch1.rxn_freqs[traj_ind, -1, rxn_ind]
                    sand.Props_integ[traj_ind, time_ind, rxn_ind] += batch1.Props_integ[traj_ind, -1, rxn_ind]
                    sand.traj_derivs[traj_ind, time_ind, rxn_ind] += batch1.traj_derivs[traj_ind, -1, rxn_ind]
                
                # Add to the end for gas phase populations
                for spec_ind in range( batch2.runAvg.Species['n_surf'] , batch2.runAvg.Species['n_surf'] + batch2.runAvg.Species['n_gas'] ):
                    sand.species_pops[traj_ind, time_ind, spec_ind] += batch1.species_pops[traj_ind, -1, spec_ind]
        
        # Combine the data
        sand.species_pops = np.concatenate( [batch1.species_pops, sand.species_pops[:,1::,:]], axis = 1 )
        sand.rxn_freqs = np.concatenate( [batch1.rxn_freqs, sand.rxn_freqs[:,1::,:]], axis = 1 )
        sand.Props_integ = np.concatenate( [batch1.Props_integ, sand.Props_integ[:,1::,:]], axis = 1 )
        sand.traj_derivs = np.concatenate( [batch1.traj_derivs, sand.traj_derivs[:,1::,:]], axis = 1 )
        
        sand.avg_updated = False
            
        return sand