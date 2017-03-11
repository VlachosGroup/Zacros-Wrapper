import os
import numpy as np
import matplotlib as mat
mat.use('Agg')
import matplotlib.pyplot as plt
import copy
import sys

from KMC_Run import KMC_Run
from Helper import FileIO, Stats
import time

'''
Has data from muliple trajectories with the same input, but different 
random seeds and possibly different initial states
'''

class Replicates:
    
    def __init__(self):
             
        # General info
        self.ParentFolder = ''
        self.runtemplate = KMC_Run()                             # Use this to build replicate jobs
        self.runAvg = KMC_Run()            # Values are averages of all runs from runList
        self.n_runs = 0
        
        # Input data for different trajectories
        self.run_dirs = []
        self.rand_seeds = []
        
        # Output data taken accross all trajectories
        # Put here various multidimensional numpy arrays which will make it easier to do statistics
        self.t_vec = []                     # times at which data is recorded
        self.species_pops = []
        self.rxn_freqs = []
        self.History_final_snaps = []
        #self.props = []            # We do not print out this output file
        self.Props_integ = []
        self.traj_derivs = []
        self.events_total = []
        self.CPU_total = []
        
        # Analysis
        self.Product = ''
        self.TOF = 0
        self.TOF_error = 0
        self.TOF_erg = 0
        self.TOF_error_erg = 0
        self.NSC_inst = []
        self.NSC_ci_inst = []
        self.NSC_erg = []
        self.NSC_ci_erg = []
    
    
    def BuildJobFiles(self, init_states = []):
        
        FileIO.ClearFolderContents(self.ParentFolder)    
        
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
                txt.write(self.runtemplate.exe_file)
            
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
                txt.write(self.runtemplate.exe_file)
    
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
                    
        sys.stdout.write('Jobs in ' + self.ParentFolder + ' have finished')
        sys.stdout.flush()
    
    def RunAllJobs_serial(self):       # Serial version of running all jobs

        # loop over array of directories and execute the Zacros executable in that folder
        for fldr in self.run_dirs:
            self.runtemplate.Path = fldr
            self.runtemplate.Run_sim()
    
    def ReadMultipleRuns(self):     # Can take ~1 minutes to use this method
        
        sys.stdout.write('Reading all runs in ' + self.ParentFolder)
        sys.stdout.flush()
        
        dummy_run = KMC_Run()       # Use this to transfer information
        
        # If directory list is empty, fill it with the directories in the current folder which have finished jobs
        if self.run_dirs == []:
        
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
            # so you will not have to store separate KMC_Run objects
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
    

    # Create a KMC run object with averaged species numbers, reaction firings, and propensities
    def AverageRuns(self):
        
        self.runAvg.Path = self.ParentFolder
        self.runAvg.Specnum['t'] = self.t_vec
        
        self.runAvg.Specnum['spec'] = np.mean(self.species_pops, axis = 0)
        self.runAvg.Procstat['events'] = np.mean(self.rxn_freqs, axis = 0)
        
        #self.runAvg.Binary['prop'] = np.mean(self.props, axis = 0)
        if not self.runAvg.Binary['propCounter'] == '':
            self.runAvg.Binary['propCounter'] = np.mean(self.Props_integ, axis = 0)
        
        self.runAvg.Performance['events_occurred'] = np.mean(self.events_total)
        self.runAvg.Performance['CPU_time'] = np.mean(self.CPU_total)
     
    def ComputeSA(self, product, window = [0.0, 1.0]):          # Need to make this compatible with irreversible reactions
        
        self.AverageRuns()      # make sure this has been done
 
        # Find index of product molecule
        try:
            product_ind = self.runAvg.Species['n_surf'] + self.runAvg.Species['gas_spec'].index(product)           # Find the index of the product species and adjust index to account for surface species
        except:
            raise Exception('Product species ' + product + ' not found.')
        
        nRxns = len(self.runAvg.Reactions['Nu'])
        TOF_stoich = np.zeros(nRxns)
        for i, elem_stoich in enumerate(self.runAvg.Reactions['Nu']):
            TOF_stoich[i] = elem_stoich[product_ind]
            
        D_stoich = np.diag(TOF_stoich)      # will use this for computing rates
        TOF_stoich.shape = [nRxns, 1]
        
        # Find indices for the beginning and end of the time window
        start_t = window[0] * self.runAvg.Specnum['t'][-1]
        end_t = window[1] * self.runAvg.Specnum['t'][-1]
        start_ind = self.runAvg.time_search(start_t)
        end_ind = self.runAvg.time_search(end_t)
        start_t = self.runAvg.Specnum['t'][start_ind]
        end_t = self.runAvg.Specnum['t'][end_t]        
        
        # Find instantaneous rates for each trajectory at each time in the window
        props_inst = ( self.Props_integ[:,end_ind,:] - self.Props_integ[:,end_ind - 1,:] ) / ( self.t_vec[end_ind] - self.t_vec[end_ind-1] )
        TOF_vec_inst = np.dot(props_inst, TOF_stoich)
        TOF_vec_inst.shape = [self.n_runs,1]
        
        props_inst_frac = np.dot(props_inst, D_stoich)
        mean_fracs = np.mean(props_inst_frac, axis = 0)
        mean_fracs = mean_fracs / np.sum(mean_fracs)            # add this to NSC
        mean_fracs = mean_fracs[::2] + mean_fracs[1::2]
        
        TOF_stats = Stats.mean_ci(TOF_vec_inst)                 # Compute the rate
        self.TOF = TOF_stats[0]
        self.TOF_error = TOF_stats[1]
        
        # Find ergodic rates for each trajectory at each time in the window
        props_erg = ( self.Props_integ[:,end_ind,:] - self.Props_integ[:,start_ind,:] )  / ( self.t_vec[end_ind] - self.t_vec[start_ind] )
        TOF_vec_erg = np.dot(props_erg, TOF_stoich)
        TOF_vec_erg.shape = [self.n_runs,1]
        
        props_erg_frac = np.dot(props_erg, D_stoich)
        mean_fracs_erg = np.mean(props_erg_frac, axis = 0)
        mean_fracs_erg = mean_fracs_erg / np.sum(mean_fracs_erg)            # add this to NSC
        mean_fracs_erg = mean_fracs_erg[::2] + mean_fracs_erg[1::2]
        
        TOF_stats_erg = Stats.mean_ci(TOF_vec_erg)          # Compute the rate
        self.TOF_erg = TOF_stats_erg[0]
        self.TOF_error_erg = TOF_stats_erg[1]
        
        # Find trajectory derivatives for each trajectory at the end of the window
        Wdata = ( self.traj_derivs[:,end_ind,::2] + self.traj_derivs[:,end_ind,1::2] )  - ( self.traj_derivs[:,start_ind,::2] + self.traj_derivs[:,start_ind,1::2] )
        
        # Perform sensitivity analysis
        W_and_rate_data = np.hstack([Wdata, TOF_vec_inst / self.TOF, TOF_vec_erg / self.TOF_erg])
        cov_out = Stats.cov_mat_ci( np.transpose(W_and_rate_data) )  
        self.NSC_inst = cov_out['cov_mat'][:-2:,-2] + mean_fracs
        self.NSC_ci_inst = cov_out['ci_mat'][:-2:,-2]
        self.NSC_erg = cov_out['cov_mat'][:-2:,-1] + mean_fracs_erg
        self.NSC_ci_erg = cov_out['ci_mat'][:-2:,-1]
        
        
    def PlotSensitivities(self): 
        
        FileIO.PlotOptions()
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        
        for i in range (self.runAvg.Reactions['nrxns']):
            cutoff = 0.05
            if self.NSC_inst[i] + self.NSC_ci_inst[i] > cutoff or self.NSC_inst[i] - self.NSC_ci_inst[i] < -cutoff:     
                plt.barh(ind-0.9, self.NSC_inst[i], width, color='r', xerr = self.NSC_ci_inst[i], ecolor='k')               
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
        
        plt.savefig(os.path.join(self.ParentFolder, 'SA_inst_output.png'))
        plt.close()
        
        plt.figure()
        width = 0.8
        ind = 0
        yvals = []
        ylabels = []
        
        for i in range (self.runAvg.Reactions['nrxns']):
            cutoff = 0.05
            if self.NSC_erg[i] + self.NSC_ci_erg[i] > cutoff or self.NSC_erg[i] - self.NSC_ci_erg[i] < -cutoff:     
                plt.barh(ind-0.9, self.NSC_erg[i], width, color='r', xerr = self.NSC_ci_erg[i], ecolor='k')               
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
        
        plt.savefig(os.path.join(self.ParentFolder, 'SA_erg_output.png'))
        plt.close()
    
    def WriteSA_output(self):
    
        with open(os.path.join(self.ParentFolder, 'SA_output.txt'), 'w') as txt:
            txt.write('Normalized sensitivity coefficients \n\n')
            txt.write('Turnover frequency: \t' + '{0:.3E} \t'.format(self.TOF) + '+- {0:.3E} \t'.format(self.TOF_error) + '\n\n')               
            txt.write('Reaction name \t NSC \t NSC confidence \n')

            for rxn_ind in range(self.runAvg.Reactions['nrxns']):
                txt.write(self.runAvg.Reactions['names'][rxn_ind] + '\t' + '{0:.3f} +- \t'.format(self.NSC_inst[rxn_ind]) + '{0:.3f}\t'.format(self.NSC_ci_inst[rxn_ind]) + '{0:.3f} +- '.format(self.NSC_erg[rxn_ind]) + '{0:.3f}'.format(self.NSC_ci_erg[rxn_ind]) + '\n')
                
    def FD_SA(self, rxn_inds = [1], pert_frac = 0.05, n_runs = 20, setup = True, exec_run = True, analyze_bool = True):     # Need to redo this function
        
        # Create objects for perturbed systems
        plus = copy.deepcopy(self)
        minus = copy.deepcopy(self)
        FD_list = [plus, minus]

        for FD in FD_list:
            FD.n_runs = n_runs
        
        # Adjust pre-exponential factors in each
        adjust_plus = np.ones(self.runtemplate.Reactions['nrxns'])
        adjust_minus = np.ones(self.runtemplate.Reactions['nrxns'])
        for rxn_ind in rxn_inds:
            adjust_plus[rxn_ind-1] = 1 + pert_frac
            adjust_minus[rxn_ind-1] = 1 / (1 + pert_frac)
        plus.runtemplate.AdjustPreExponentials(adjust_plus)
        minus.runtemplate.AdjustPreExponentials(adjust_minus)
        
        # Set subfolder for perturbed runs
        plus.ParentFolder = self.ParentFolder + 'plus'
        minus.ParentFolder = self.ParentFolder + 'minus'
        
        if setup:
        
            for FD in FD_list:
                for rxn_ind in rxn_inds:
                    FD.ParentFolder = FD.ParentFolder + '_' + str(rxn_ind)
                FD.ParentFolder = FD.ParentFolder + '/'
            
                # Build folders for runs
                if not os.path.exists(FD.ParentFolder):
                        os.makedirs(FD.ParentFolder)
                FD.BuildJobs()
        
        if exec_run:
            for FD in FD_list:
                FD.RunAllJobs()
        
        ''' Analyze results '''
        if analyze_bool:
            for FD in FD_list:
                FD.ReadMultipleRuns()
                FD.AverageRuns()
                
            plus_stats = []
            for run in plus.runList:
                plus_stats.append(run.ComputeTOF(self.Product)['TOF'])
                
            minus_stats = []
            for run in minus.runList:
                minus_stats.append(run.ComputeTOF(self.Product)['TOF'])
    
            all_TOFs = plus_stats + minus_stats
            TOF_mean = np.mean(all_TOFs)
            diff_stats = Stats.diff_ci(plus_stats, minus_stats)   
            
            NSC = diff_stats[0] / TOF_mean / (2 * pert_frac)
            NSC_ci = diff_stats[1] / TOF_mean / (2 * pert_frac)
            
            return [NSC, NSC_ci]
            
            
    @staticmethod
    def time_sandwich(batch1, batch2):
        
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
            
        return sand