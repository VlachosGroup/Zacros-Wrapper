import os
import numpy as np
import matplotlib as mat
#mat.use('Agg')
import matplotlib.pyplot as plt
import copy
import sys

from KMC_Run import kmc_traj
from Helper import *
import time



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
        self.N_batches = 1000
        self.Nbpt = None
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
        
        
    def Compute_ACF(self):          # Need to also put error bars on this
    
        '''
        Check steady state on the basis of batch means
        '''
        
        if not self.avg_updated:
            self.AverageRuns()
        
        # Find index of product molecule
        try:
            product_ind = self.runAvg.Species['n_surf'] + self.runAvg.Species['gas_spec'].index(self.gas_product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species ' + self.gas_product + ' not found.')
            
        nRxns = len(self.runAvg.Reactions['Nu'])
        TOF_stoich = np.zeros(nRxns)
        for i, elem_stoich in enumerate(self.runAvg.Reactions['Nu']):
            TOF_stoich[i] = elem_stoich[product_ind]
            
        self.Nbpt = np.max( [ 3 , (self.N_batches-1) / self.n_runs + 1 ] )            # Set the number of batches per trajectory
        
        bin_edges = np.linspace(0, self.t_vec[-1], self.Nbpt + 1)
        rate_data = np.zeros([self.n_runs, self.Nbpt])
        
        for i in range(self.Nbpt):
        
            idb_start = self.runAvg.time_search_interp(bin_edges[i])
            idb_end = self.runAvg.time_search_interp(bin_edges[i+1])
            
            prop_integ_start = idb_start[1][0] * self.Props_integ[:, idb_start[0][0], :] + idb_start[1][1] * self.Props_integ[:, idb_start[0][1], :]
            prop_integ_end = idb_end[1][0] * self.Props_integ[:, idb_end[0][0], :] + idb_end[1][1] * self.Props_integ[:, idb_end[0][1], :]
            
            rate_data[:,i] = np.dot ( ( prop_integ_end - prop_integ_start ) / ( bin_edges[i+1] - bin_edges[i] ) , TOF_stoich )
        
        # Compute the autocorrelation function for the batch means of the rate
        c1 = rate_data[:,1:-1:]
        c1 = c1.reshape( np.prod(c1.shape) )
        c2 = rate_data[:,2::]
        c2 = c2.reshape( np.prod(c2.shape) )
        
        if np.var(c1) * np.var(c2) == 0:        # Covers the case where we have zero rate
            return 1.0
        else:
            c = ( np.mean(c1 * c2) - np.mean(c1) * np.mean(c2) ) / np.sqrt( np.var(c1) * np.var(c2) )
            return c
        
    
    def PerformSA(self, delta_t = None, ergodic = False,):          # Need to make this compatible with irreversible reactions
        
        '''
        Perform sensitivity analysis with a combination of
        time and trajectory averaging
        
        delta_t             size of the time window
        ergodic             True: use the rate at the end of the time interval, False: average the rate over the entire time interval 
        
        The time-averaging could be improved in here. Right now it is a little weird...
        You can use linear interpolation to get data between time points
        Can also make it so that the statistics make sense for non-evenly spaced time points
        '''
        
        # Use entire trajectory length as the default time window
        if delta_t is None:
            delta_t = self.t_vec[-1] / self.Nbpt
        
        if not self.avg_updated:
            self.AverageRuns()
 
        # Find index of product molecule
        try:
            product_ind = self.runAvg.Species['n_surf'] + self.runAvg.Species['gas_spec'].index(self.gas_product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species ' + self.gas_product + ' not found.')
        
        nRxns = len(self.runAvg.Reactions['Nu'])
        TOF_stoich = np.zeros(nRxns)
        for i, elem_stoich in enumerate(self.runAvg.Reactions['Nu']):
            TOF_stoich[i] = elem_stoich[product_ind]
            
        if delta_t > self.t_vec[-1] or delta_t < 0:
            print 'Time window: ' + str(delta_t)
            print 'Final time: ' + str(self.t_vec[-1])
            raise Exception('Time window is too large. Insufficient sampling.')
        
        delt_ind = self.runAvg.time_search(delta_t)
        n_data_points = ( len( self.t_vec ) - delt_ind ) * self.n_runs
        
        NSCs = np.zeros(self.traj_derivs.shape[2])
        W_mean = np.zeros(self.traj_derivs.shape[2])
        mean_rate = 0
        rate_contributions_avgs = np.zeros(self.traj_derivs.shape[2])
        
        for traj_ind in range(self.n_runs):
            print 'Trajectory #: ' + str(traj_ind+1)
            for time_ind in range(delt_ind, len( self.t_vec ) ):
            
                t = self.t_vec[time_ind]
                t_begin_win_ind = self.runAvg.time_search(t - delta_t)
                t_begin_win = self.t_vec[t_begin_win_ind]
                
                if ergodic:             # Average the rate over the entire time window
                    rate_contributions = ( self.Props_integ[traj_ind, time_ind,:] - self.Props_integ[traj_ind, t_begin_win_ind, :] ) / ( self.t_vec[time_ind] - self.t_vec[t_begin_win_ind] ) * TOF_stoich
                else:                   # Take the final value of the rate
                    rate_contributions = ( self.Props_integ[traj_ind, time_ind,:] - self.Props_integ[traj_ind, time_ind-1, :] ) / ( self.t_vec[time_ind] - self.t_vec[time_ind-1] ) * TOF_stoich 
                
                rate_contributions_avgs = rate_contributions_avgs + rate_contributions
                rate = np.sum(rate_contributions)
                
                Ws = self.traj_derivs[traj_ind, time_ind , :] - self.traj_derivs[ traj_ind, t_begin_win_ind , :]
                
                mean_rate += rate
                W_mean = W_mean + Ws
                NSCs = NSCs + Ws * rate
        
        # Finalize the averaging
        NSCs = NSCs / n_data_points
        mean_rate = mean_rate / n_data_points
        W_mean = W_mean / n_data_points
        rate_contributions_avgs = rate_contributions_avgs / n_data_points
        
        NSCs = ( NSCs - W_mean * mean_rate ) / mean_rate + rate_contributions_avgs / mean_rate
        
        # Combine forward and backwards NSCs for reversible reactions
        NSC_rev = np.zeros(self.runAvg.Reactions['nrxns'])
        ind = 0
        for i in range(self.runAvg.Reactions['nrxns']):
            
            if self.runAvg.Reactions['is_reversible'][i]:
                NSC_rev[i] = NSCs[ind] + NSCs[ind+1]
                ind += 2
            else:
                NSC_rev[i] = NSCs[ind]
                ind += 1
        
        print 'Normalized sensitivity coefficients: ' + str(NSC_rev)
        print 'Turnover frequency: ' + str(mean_rate)
        
        self.NSC = NSC_rev
        self.NSC_ci = np.zeros([len(self.NSC)])
        
        
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