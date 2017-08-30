import numpy as np
import os
import copy

from Replicates import Replicates
from Helper import *

import matplotlib as mat
import matplotlib.pyplot as plt
  
def ReachSteadyStateAndRescale(kmc_template, scale_parent_fldr, n_runs = 16, n_batches = 1000, 
                                acf_cut = 0.05, include_stiff_reduc = True, max_events = int(1e3), 
                                max_iterations = 30, ss_inc = 2.0, n_samples = 100, platform = 'Squidward'):

    '''
    Handles rate rescaling and continuation of KMC runs
    
    kmc_template           :     kmc_traj object with information about the physical system
    scale_parent_fldr      :     Working directory
    n_runs                 :     Number of trajectories to run, also the number of processors
    include_stiff_reduc    :     True to allow for scaledown, False to turn this feature off
    max_events             :     Maximum number of events for the first iteration
    max_iterations         :     Maximum number of iterations to use
    ss_inc                 :     Factor to scale the final time by if you have not yet reached steady state
    n_samples              :     Number of time points to sample for each trajectory
    platform               :     Squidward or Farber
    '''
    
    prev_batch = Replicates()       # Set this if the starting iteration is not 1
    initial_states = []
    
    # Placeholder variables
    if not os.path.exists(scale_parent_fldr):
        os.makedirs(scale_parent_fldr)
    ClearFolderContents(scale_parent_fldr)
    SDF_vec = None        # scaledown factors for each iteration
    
    # Convergence variables
    is_steady_state = False
    iteration = 1              
    
    scale_final_time = ss_inc

    while not is_steady_state and iteration <= max_iterations:
        
        # Make folder for iteration
        iter_fldr = os.path.join(scale_parent_fldr, 'Iteration_' + str(iteration))
        if not os.path.exists(iter_fldr):
            os.makedirs(iter_fldr)
            
        # Create object for batch
        cur_batch = Replicates()
        cur_batch.ParentFolder = iter_fldr
        cur_batch.n_runs = n_runs
        cur_batch.N_batches = n_batches
        cur_batch.Set_kmc_template(kmc_template)        # Set template KMC trajectory
        
        if iteration == 1:              # Sample on events, because we do not know the time scales
        
            # Set sampling parameters
            cur_batch.runtemplate.simin.MaxStep = max_events
            cur_batch.runtemplate.simin.SimTime_Max = 'inf'
            cur_batch.runtemplate.simin.WallTime_Max = 'inf'
            cur_batch.runtemplate.simin.restart = False
            
            cur_batch.runtemplate.simin.procstat = ['event', np.max( [max_events / n_samples, 1] ) ]
            cur_batch.runtemplate.simin.specnum = ['event', np.max( [max_events / n_samples, 1] ) ]
            cur_batch.runtemplate.simin.hist = ['event', np.max( [max_events * (n_samples-1) / n_samples, 1] )]       # only record the initial and final states

            SDF_vec = np.ones( cur_batch.runtemplate.mechin.get_num_rxns() )         # Initialize scaledown factors
        
        elif iteration > 1:             # Time sampling

            # Change sampling
            cur_batch.runtemplate.simin.MaxStep = 'inf'
            cur_batch.runtemplate.simin.WallTime_Max = 'inf'
            cur_batch.runtemplate.simin.restart = False
            cur_batch.runtemplate.simin.SimTime_Max = prev_batch.t_vec[-1] * scale_final_time
            cur_batch.runtemplate.simin.SimTime_Max = float('{0:.3E} \t'.format( cur_batch.runtemplate.simin.SimTime_Max ))     # round to 4 significant figures
            cur_batch.runtemplate.simin.procstat = ['time', cur_batch.runtemplate.simin.SimTime_Max / n_samples]
            cur_batch.runtemplate.simin.specnum = ['time', cur_batch.runtemplate.simin.SimTime_Max / n_samples]
            cur_batch.runtemplate.simin.hist = ['time', cur_batch.runtemplate.simin.SimTime_Max ]
            
            # Adjust pre-exponential factors based on the stiffness assessment of the previous iteration
            if include_stiff_reduc:
                cur_batch.runtemplate.AdjustPreExponentials(SDF_vec)
            
            # Use continuation
            initial_states = prev_batch.History_final_snaps
        
        # Run jobs and read output
        cur_batch.BuildJobFiles(init_states = initial_states)
        cur_batch.RunAllJobs_parallel_JobArray(server = platform, job_name = 'Iteration_' + str(iteration) )
        cur_batch.ReadMultipleRuns()
        
        if iteration == 1:
            cum_batch = copy.deepcopy(cur_batch)
        else:
            cum_batch = Replicates.time_sandwich(prev_batch, cur_batch)         # combine with previous data          
        
        # Test steady-state
        cum_batch.AverageRuns()
        acf_data = cum_batch.Compute_ACF()
        is_steady_state = False
        if not acf_data['ACF'] is None:
            if acf_data['ACF'] + acf_data['ACF_ci'] < 0.08:
                is_steady_state = True
        
        # Make sure the error in the rate is less than 5% of the rate
        if acf_data['Rate'] == 0:
            is_steady_state = False
        else:
            if acf_data['Rate_ci'] / acf_data['Rate'] > 0.05:
                is_steady_state = False
        
        # Record information about the iteration
        cum_batch.runAvg.PlotGasSpecVsTime()
        cum_batch.runAvg.PlotSurfSpecVsTime()
        
        cur_batch.AverageRuns()
        cur_batch.runAvg.PlotElemStepFreqs()
        scaledown_data = ProcessStepFreqs(cur_batch.runAvg)         # compute change in scaledown factors based on simulation result
        delta_sdf = scaledown_data['delta_sdf']

        ## Record iteartion data in output file
        #with open(os.path.join(iter_fldr, 'Iteration_' + str(iteration) + '_summary.txt'), 'w') as txt:  
        #    
        #    txt.write('----- Iteration #' + str(iteration) + ' -----\n')
        #    txt.write('t_final: {0:.3E} \n'.format(cum_batch.runAvg.Specnum['t'][-1]))
        #    txt.write('steady-state: ' + str(is_steady_state) + '\n')
        #    
        #    txt.write('Reaction name \t')
        #    txt.write('Scaledown factor \t')
        #    txt.write('Reaction speed \t')
        #    txt.write('Total firings \t')
        #    txt.write('Net firings \n')
        #
        #    for rxn_ind in range(len(cum_batch.runAvg.Reactions['names'])):
        #        txt.write(cum_batch.runAvg.Reactions['names'][rxn_ind] + '\t')
        #        txt.write('{0:.3E} \t'.format(delta_sdf[rxn_ind]))
        #        txt.write(scaledown_data['rxn_speeds'][rxn_ind] + '\t')
        #        txt.write('{0:.3E} \t'.format(scaledown_data['tot'][rxn_ind]))
        #        txt.write('{0:.3E} \n'.format(scaledown_data['net'][rxn_ind]))
        
        # Update scaledown factors
        for ind in range(len(SDF_vec)):
            SDF_vec[ind] = SDF_vec[ind] * delta_sdf[ind]
            
        scale_final_time = np.max( [1.0/np.min(delta_sdf), ss_inc] )
        
        prev_batch = copy.deepcopy(cum_batch)
        iteration += 1
        #is_steady_state = False
    return cum_batch
    

def ProcessStepFreqs(run, stiff_cut = 100.0, equilib_cut = 0.1):        # Change to allow for irreversible reactions
    
    '''
    Takes an average KMC trajectory and assesses the reaction frequencies to identify fast reactions
    Process KMC output and determine how to further scale down reactions
    Uses algorithm from A. Chatterjee, A.F. Voter, Accurate acceleration of kinetic Monte Carlo simulations through the modification of rate constants, J. Chem. Phys. 132 (2010) 194101.
    '''
    
    delta_sdf = np.ones( run.mechin.get_num_rxns() )    # initialize the marginal scaledown factors
    rxn_speeds = []
    
    # data analysis
    freqs = run.procstatout.events[-1,:]
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
        #alpha_UB = N_f .* delta ./ log(1 ./ delta) + 1             # Chatterjee formula
        delta_sdf[i] = np.min([1.0, stiff_cut / N_f])
        
    return {'delta_sdf': delta_sdf, 'rxn_speeds': rxn_speeds, 'tot': tot_freqs, 'net': net_freqs}
    
    
def ReadScaledown(RunPath, plot_analysis = False, fldrs_cut = None, verbose = False, product = None, n_batches = 1000):
    
    '''
    Read a scaledown that has already been run
    '''

    # Prepare data to graph
    batch_lengths = []
    rates_vec = []
    rates_vec_ci = []
    acf_vec = []
    acf_vec_ci = []
    
    
    # Count the iterations
    n_folders = len(os.listdir(RunPath))
    if not fldrs_cut is None:
       
        n_folders = fldrs_cut

    print str(n_folders) + ' iterations found'

    cum_batch = None
    for ind in range(1,n_folders+1):
        
        x = Replicates()
        x.ParentFolder = os.path.join(RunPath, 'Iteration_' + str(ind))
        x.ReadMultipleRuns()

        if ind == 1:
            cum_batch = x
        else:
            cum_batch = Replicates.time_sandwich(cum_batch, x)
            
        cum_batch.N_batches = n_batches
        cum_batch.gas_product = product
        acf_data = cum_batch.Compute_ACF()
        print 'Iteration ' + str(ind)
        print 'Batches per trajectory: ' + str(cum_batch.Nbpt)
        print 'Batch length ' + str(cum_batch.batch_length)
        print 'Rate: ' + str(acf_data['Rate'])
        print 'Rate CI: ' + str(acf_data['Rate_ci'])
        print 'ACF: ' + str(acf_data['ACF'])
        print 'ACF CI: ' + str(acf_data['ACF_ci'])

        print '\n'
        
        batch_lengths.append(cum_batch.batch_length)
        rates_vec.append(acf_data['Rate'])
        rates_vec_ci.append(acf_data['Rate_ci'])
        acf_vec.append(acf_data['ACF'])
        acf_vec_ci.append(acf_data['ACF_ci'])
            
    print batch_lengths
    print rates_vec
    print rates_vec_ci
    print acf_vec
    print acf_vec_ci
        
            
    return cum_batch

    