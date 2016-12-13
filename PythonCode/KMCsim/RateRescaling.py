# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 20:28:48 2016

@author: RDX
"""

import numpy as np
import os
import copy

from Replicates import Replicates
from KMC_Run import KMC_Run
from Helper import Helper

class RateRescaling:
    
    def __init__(self):
        
        self.summary_filename = 'rescaling_output.txt'
        self.scale_parent_fldr = ''
        self.batch = Replicates()
        
    def ReachSteadyStateAndRescale(self, Product, template_folder, exe, include_stiff_reduc = True, max_events = int(1e4), max_iterations = 15, stiff_cutoff = 1, ss_inc = 2.0, n_samples = 100, n_runs = 10):

        Helper.ClearFolderContents(self.scale_parent_fldr)

        # Placeholder variables
        prev_batch = Replicates()
        
        # Convergence variables
        is_steady_state = False
        unstiff = False
        converged = unstiff and is_steady_state
        iteration = 1
        SDF_vec = []        # scaledown factors for each iteration
        scale_final_time = ss_inc

        while not converged and iteration <= max_iterations:
            
            # Make folder for iteration
            iter_fldr = self.scale_parent_fldr + 'Iteration_' + str(iteration) + '/'
            if not os.path.exists(iter_fldr):
                os.makedirs(iter_fldr)
                
            # Create object for batch
            cur_batch = Replicates()
            cur_batch.ParentFolder = iter_fldr
            cur_batch.n_runs = n_runs
            cur_batch.Product = Product
            
            cur_batch.runtemplate = KMC_Run()
            cur_batch.runtemplate.exe_file = exe
            cur_batch.runtemplate.Path = template_folder
            cur_batch.runtemplate.ReadAllInput()
            
            if iteration == 1:              # Event sampling
            
                # Set sampling parameters
                cur_batch.runtemplate.Conditions['MaxStep'] = max_events
                cur_batch.runtemplate.Conditions['SimTime']['Max'] = 'inf'
                cur_batch.runtemplate.Conditions['WallTime']['Max'] = 'inf'
                cur_batch.runtemplate.Conditions['restart'] = False
                
                cur_batch.runtemplate.Report['procstat'] = ['event', max_events / n_samples]
                cur_batch.runtemplate.Report['specnum'] = ['event', max_events / n_samples]
                cur_batch.runtemplate.Report['hist'] = ['event', max_events * (n_samples-1) / n_samples]       # only record the initial and final states

                SDF_vec = np.ones(cur_batch.runtemplate.Reactions['nrxns'])         # Initialize scaledown factors
            
            elif iteration > 1:             # Time sampling

                # Change sampling
                cur_batch.runtemplate.Conditions['MaxStep'] = 'inf'
                cur_batch.runtemplate.Conditions['WallTime']['Max'] = 'inf'
                cur_batch.runtemplate.Conditions['restart'] = False
                cur_batch.runtemplate.Conditions['SimTime']['Max'] = prev_batch.runList[0].Performance['t_final'] * scale_final_time
                cur_batch.runtemplate.Conditions['SimTime']['Max'] = float('{0:.3E} \t'.format(cur_batch.runtemplate.Conditions['SimTime']['Max']))     # round to 4 significant figures
                cur_batch.runtemplate.Report['procstat'] = ['time', cur_batch.runtemplate.Conditions['SimTime']['Max'] / n_samples]
                cur_batch.runtemplate.Report['specnum'] = ['time', cur_batch.runtemplate.Conditions['SimTime']['Max'] / n_samples]
                cur_batch.runtemplate.Report['hist'] = ['time', cur_batch.runtemplate.Conditions['SimTime']['Max']]
                
                if include_stiff_reduc:
                    cur_batch.runtemplate.AdjustPreExponentials(SDF_vec)
                
            cur_batch.BuildJobsFromTemplate()
                
            # Use continuation
            if iteration > 1:
                for run_ind in range(n_runs):
                    cur_batch.runList[run_ind].StateInput['Type'] = 'history'
                    cur_batch.runList[run_ind].StateInput['Struct'] = prev_batch.runList[run_ind].History[-1]
            
            # Run jobs and read output
            cur_batch.BuildJobFiles()
            cur_batch.SubmitJobArray()
            cur_batch.WaitForJobs()
            cur_batch.ReadMultipleRuns()
            
            # Test steady-state
            if iteration == 1:
                is_steady_state = False
            else:
                cur_batch.AverageRuns()
                cur_batch.ParentFolder = iter_fldr
                cur_batch.runAvg.Path = iter_fldr
                correl = cur_batch.CheckAutocorrelation(Product)
                is_steady_state = np.abs(correl[0]) < 0.05
                
                # Record information about the iteration
                cur_batch.runAvg.CalcRateTraj(Product)
        
                cur_batch.runAvg.PlotSurfSpecVsTime()        
                cur_batch.runAvg.PlotIntPropsVsTime()
                cur_batch.runAvg.PlotRateVsTime()  
            
            # Test stiffness
            cur_batch.AverageRuns()
            scaledown_data = RateRescaling.ProcessStepFreqs(cur_batch.runAvg)         # compute change in scaledown factors based on simulation result
            delta_sdf = scaledown_data['delta_sdf']
            rxn_speeds = scaledown_data['rxn_speeds']
            if include_stiff_reduc:
                unstiff = np.max(np.abs(np.log10(delta_sdf))) < stiff_cutoff
            else:
                unstiff = True
    
            # Record iteartion data in output file
            with open(iter_fldr + 'Iteration_summary.txt', 'w') as txt:   
                txt.write('----- Iteration #' + str(iteration) + ' -----\n')
                txt.write('t_final: {0:.3E} \n'.format(cur_batch.runAvg.Specnum['t'][-1]))
                txt.write('stiff: ' + str(not unstiff) + '\n')
                txt.write('steady-state: ' + str(is_steady_state) + '\n')
                for rxn_name in cur_batch.runAvg.Reactions['names']:
                    txt.write(rxn_name + '\t')
                txt.write('\n')
                for sdf in delta_sdf:
                    txt.write('{0:.3E} \t'.format(sdf))
                txt.write('\n')
                for rxn_speed in rxn_speeds:
                    txt.write(rxn_speed + '\t')
            
            # Update scaledown factors
            for ind in range(len(SDF_vec)):
                SDF_vec[ind] = SDF_vec[ind] * delta_sdf[ind]
                
            scale_final_time = np.max( [1.0/np.min(delta_sdf), ss_inc] )
            
            prev_batch = copy.deepcopy(cur_batch)
            converged = unstiff and is_steady_state
            iteration += 1
            
        cur_batch.ComputeStats(Product)
        cur_batch.WriteSA_output()
        cur_batch.PlotSensitivities()
    
    # Process KMC output and determine how to further scale down reactions
    @staticmethod
    def ProcessStepFreqs(run, stiff_cut = 100, equilib_cut = 0.05):
        
        delta_sdf = np.ones(run.Reactions['nrxns'])    # initialize the marginal scaledown factors
        rxn_speeds = []
        
        # data analysis
        freqs = run.Procstat['events'][-1,:]
        fwd_freqs = freqs[0::2]
        bwd_freqs = freqs[1::2]
        net_freqs = fwd_freqs - bwd_freqs
        tot_freqs = fwd_freqs + bwd_freqs
        
        fast_rxns = []
        slow_rxns = []        
        for i in range(len(tot_freqs)):
            if tot_freqs[i] == 0:
                slow_rxns.append(i)
                rxn_speeds.append('slow')
            else:
                PE = float(net_freqs[i]) / tot_freqs[i]
                if np.abs(PE) < equilib_cut:
                    fast_rxns.append(i)
                    rxn_speeds.append('fast')
                else:
                    slow_rxns.append(i)
                    rxn_speeds.append('slow')
        
        # Find slow scale rate
        slow_freqs = [1.0]      # put an extra 1 in case no slow reactions occur
        for i in slow_rxns:
            slow_freqs.append(tot_freqs[i])
        slow_scale = np.max(slow_freqs)
        
        # Adjust fast reactions closer to the slow scale
        for i in fast_rxns:
            delta_sdf[i] = np.min([1.0, stiff_cut * float(slow_scale) / tot_freqs[i]])
            
        return {'delta_sdf': delta_sdf, 'rxn_speeds': rxn_speeds}