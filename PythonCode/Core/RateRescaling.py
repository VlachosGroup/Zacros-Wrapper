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
from Helper import FileIO

class RateRescaling:
    
    def __init__(self):
        
        self.summary_filename = 'rescaling_output.txt'
        self.scale_parent_fldr = ''
        self.batch = Replicates()
        
    def ReachSteadyStateAndRescale(self, Product, gas_stoich, template_folder, exe, include_stiff_reduc = True, max_events = int(1e4), max_iterations = 15, stiff_cutoff = 1, ss_inc = 2.0, n_samples = 100, n_runs = 10, start_iter = 1, platform = 'Farber'):

        prev_batch = Replicates()       # Set this if the starting iteration is not 1
        initial_states = []
        
        # Placeholder variables
        if start_iter == 1:
            FileIO.ClearFolderContents(self.scale_parent_fldr)
            SDF_vec = []        # scaledown factors for each iteration
        else:
            for i in range(start_iter-1):
                batch_i = Replicates()
                batch_i.Path = os.path.join(self.scale_parent_fldr, 'Iteration_' + str(i+1))
                batch_i.ReadMultipleRuns()
                prev_batch = Replicates.time_sandwich(prev_batch, batch_i)
                
                if i == start_iter-2:
                    batch_i.AverageRuns()
                    scaledown_data = RateRescaling.ProcessStepFreqs(batch_i.runAvg)         # compute change in scaledown factors based on simulation result
                    delta_sdf = scaledown_data['delta_sdf']
                    rxn_speeds = scaledown_data['rxn_speeds']                
                
            SDF_vec = prev_batch.runList[0].scaledown_factors
            # Update scaledown factors
            for ind in range(len(SDF_vec)):
                SDF_vec[ind] = SDF_vec[ind] * delta_sdf[ind]
        
        # Convergence variables
        is_steady_state = False
        unstiff = False
        converged = unstiff and is_steady_state
        iteration = start_iter              # default is 1, but can be higher and use previous data
        
        scale_final_time = ss_inc

        while not converged and iteration <= max_iterations:
            
            # Make folder for iteration
            iter_fldr = os.path.join(self.scale_parent_fldr, 'Iteration_' + str(iteration))
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
                
                # Use continuation
                initial_states = prev_batch.History_final_snaps
            
            # Run jobs and read output
            cur_batch.BuildJobFiles(init_states = initial_states)
            cur_batch.SubmitJobArray(server = platform)
            cur_batch.ReadMultipleRuns()
            cum_batch = Replicates.time_sandwich(prev_batch, cur_batch)         # combine with previous data          
            
            # Test steady-state
            cum_batch.AverageRuns()
            cum_batch.runAvg.gas_stoich = gas_stoich        # stoichiometry of gas-phase reaction
            cum_batch.runAvg.calc_net_rxn()            
            net_rxn_converged = cum_batch.runAvg.CheckNetRxnConvergence()
            product_converged = cum_batch.runAvg.CheckProductConvergence(Product)
            is_steady_state = net_rxn_converged and product_converged
            
            # Record information about the iteration
            cum_batch.runAvg.PlotGasSpecVsTime()
            cum_batch.runAvg.PlotNetGasRxnVsTime()
            cum_batch.runAvg.PlotSurfSpecVsTime()
            
            # Test stiffness
            cur_batch.AverageRuns()
            cur_batch.runAvg.PlotElemStepFreqs()
            scaledown_data = RateRescaling.ProcessStepFreqs(cur_batch.runAvg)         # compute change in scaledown factors based on simulation result
            delta_sdf = scaledown_data['delta_sdf']
            if include_stiff_reduc:
                unstiff = np.max(np.abs(np.log10(delta_sdf))) < stiff_cutoff
            else:
                unstiff = True
    
            # Record iteartion data in output file
            with open(os.path.join(iter_fldr, 'Iteration_' + str(iteration) + '_summary.txt'), 'w') as txt:  
                
                txt.write('----- Iteration #' + str(iteration) + ' -----\n')
                txt.write('t_final: {0:.3E} \n'.format(cum_batch.runAvg.Specnum['t'][-1]))
                txt.write('stiff: ' + str(not unstiff) + '\n')
                txt.write('steady-state: ' + str(is_steady_state) + '\n')
                
                txt.write('Reaction name \t')
                txt.write('Scaledown factor \t')
                txt.write('Reaction speed \t')
                txt.write('Total firings \t')
                txt.write('Net firings \n')
                
                for rxn_ind in range(len(cum_batch.runAvg.Reactions['names'])):
                    txt.write(cum_batch.runAvg.Reactions['names'][rxn_ind] + '\t')
                    txt.write('{0:.3E} \t'.format(delta_sdf[rxn_ind]))
                    txt.write(scaledown_data['rxn_speeds'][rxn_ind] + '\t')
                    txt.write('{0:.3E} \t'.format(scaledown_data['tot'][rxn_ind]))
                    txt.write('{0:.3E} \n'.format(scaledown_data['net'][rxn_ind]))
            
            # Update scaledown factors
            for ind in range(len(SDF_vec)):
                SDF_vec[ind] = SDF_vec[ind] * delta_sdf[ind]
                
            scale_final_time = np.max( [1.0/np.min(delta_sdf), ss_inc] )
            
            prev_batch = copy.deepcopy(cum_batch)
            converged = unstiff and is_steady_state
            iteration += 1
            
        return cur_batch
    
    # Process KMC output and determine how to further scale down reactions
    # Uses algorithm from A. Chatterjee, A.F. Voter, Accurate acceleration of kinetic Monte Carlo simulations through the modification of rate constants, J. Chem. Phys. 132 (2010) 194101.
    @staticmethod
    def ProcessStepFreqs(run, stiff_cut = 40.0, equilib_cut = 0.05):
        
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
            N_f = tot_freqs[i] / float(slow_scale)              # number of fast events per rare event
#            alpha_UB = N_f .* delta ./ log(1 ./ delta) + 1             # Chatterjee formula
            delta_sdf[i] = np.min([1.0, stiff_cut / N_f])
            
        return {'delta_sdf': delta_sdf, 'rxn_speeds': rxn_speeds, 'tot': tot_freqs, 'net': net_freqs}