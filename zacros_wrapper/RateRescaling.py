import numpy as np
import os
import copy

from zacros_wrapper.Replicates import *
from zacros_wrapper.utils import *

import matplotlib as mat
import matplotlib.pyplot as plt

try:
    from mpi4py import MPI      # Import MPI if available
except:
    pass


def ReachSteadyStateAndRescale(kmc_template, scale_parent_fldr, n_runs = 16, n_batches = 1000,
                                acf_cut = 0.05, include_stiff_reduc = True, max_events = int(1e3),
                                max_iterations = 30, ss_inc = 1.0, n_samples = 100, parallel_mode = 'Squidward',
                                ACF_tol = 0.08, rate_tol = 0.05):

    '''
    Handles rate rescaling and continuation of KMC runs

    :param kmc_template:            kmc_traj object with information about the physical system

    :param scale_parent_fldr:       Working directory

    :param n_runs:                  Number of trajectories to run, also the number of processors

    :param include_stiff_reduc:     True to allow for scaledown, False to turn this feature off

    :param max_events:              Maximum number of events for the first iteration

    :param max_iterations:          Maximum number of iterations to use

    :param ss_inc:                  Factor to scale the final time by if you have not yet reached steady state

    :param n_samples:               Number of time points to sample for each trajectory

    :param parallel_mode:           MPI, Squidward, or Farber
    '''

    if parallel_mode == 'MPI':
        try:
            COMM = MPI.COMM_WORLD
            COMM.Barrier()
        except:
            raise NameError('mpi4py dependency has not been imported.')

    prev_batch = Replicates()       # Set this if the starting iteration is not 1
    initial_states = None

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
        cur_batch.n_trajectories = n_runs
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
        if parallel_mode == 'MPI':
            if COMM.rank == 0:
                cur_batch.BuildJobFiles(init_states = initial_states)

            # Collect whatever has to be done in a list. Here we'll just collect a list of
            # numbers. Only the first rank has to do this.
            if COMM.rank == 0:
                jobs = cur_batch.run_dirs
                jobs = [jobs[_i::COMM.size] for _i in range(COMM.size)]             # Split into however many cores are available.
            else:
                jobs = None

            jobs = COMM.scatter(jobs, root=0)           # Scatter jobs across cores.

            # Now each rank just does its jobs and collects everything in a results list.
            # Make sure to not use super big objects in there as they will be pickled to be
            # exchanged over MPI.
            for job in jobs:
                cur_batch.runtemplate.Path = job
                cur_batch.runtemplate.Run_sim()

        else:
            cur_batch.BuildJobFiles(init_states = initial_states)
            cur_batch.RunAllTrajectories_JobArray(server = parallel_mode, job_name = 'Iteration_' + str(iteration) )

        cur_batch.ReadMultipleRuns()

        if iteration == 1:
            cum_batch = copy.deepcopy(cur_batch)
        else:
            cum_batch = append_replicates(prev_batch, cur_batch)         # combine with previous data

        # Test steady-state
        cum_batch.AverageRuns()
        acf_data = cum_batch.Compute_rate()

        print('\nIteration ' + str(iteration))
        print('Batches per trajectory: ' + str(cum_batch.Nbpt))
        print('Batch length (s): ' + str(cum_batch.batch_length))
        print('Rate: ' + str(cum_batch.rate))
        print('Rate confidence interval: ' + str(cum_batch.rate_CI))
        print('Autocorrelation: ' + str(cum_batch.ACF))
        print('Autocorrelation confidence: ' + str(cum_batch.ACF_CI))

        # Test if autocorrelation function has converged

        if cum_batch.ACF is None:
            decorrelated = False
        else:
            decorrelated = ( cum_batch.ACF + cum_batch.ACF_CI < ACF_tol)

        # Test if rate is computed with sufficient accuracy
        if cum_batch.rate == 0:
            rate_accurate = False
        else:
            rate_accurate = (cum_batch.rate_CI / cum_batch.rate < rate_tol)

        print('Decorrelated? ' + str(decorrelated))
        print('Rate accurate? ' + str(rate_accurate))
        print('\n')

        is_steady_state = decorrelated and rate_accurate

        # Record information about the iteration
        cum_batch.runAvg.PlotGasSpecVsTime()
        cum_batch.runAvg.PlotSurfSpecVsTime()

        cur_batch.AverageRuns()
        cur_batch.runAvg.PlotElemStepFreqs()
        scaledown_data = ProcessStepFreqs(cur_batch.runAvg)         # compute change in scaledown factors based on simulation result
        delta_sdf = scaledown_data['delta_sdf']

        # Update scaledown factors
        for ind in range(len(SDF_vec)):
            SDF_vec[ind] = SDF_vec[ind] * delta_sdf[ind]

        scale_final_time = np.max( [1.0/np.min(delta_sdf), ss_inc] )

        prev_batch = copy.deepcopy(cum_batch)
        iteration += 1

    return cum_batch


def ProcessStepFreqs(run, stiff_cut = 100.0, delta = 0.05, equilib_cut = 0.1):        # Change to allow for irreversible reactions

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
        #alpha_UB = N_f * delta / np.log(1 / delta) + 1             # Chatterjee formula

        #delta_sdf[i] = np.min([1.0, np.max([stiff_cut / N_f, 1. / alpha_UB ]) ])
        delta_sdf[i] = np.min([1.0, stiff_cut / N_f ])

    return {'delta_sdf': delta_sdf, 'rxn_speeds': rxn_speeds, 'tot': tot_freqs, 'net': net_freqs}


def ReadScaledown(RunPath, fldrs_cut = None, product = None, n_batches = 1000):

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

    print(str(n_folders) + ' iterations found')

    cum_batch = None
    for ind in range(1,n_folders+1):

        x = Replicates()
        x.ParentFolder = os.path.join(RunPath, 'Iteration_' + str(ind))
        x.ReadMultipleRuns()

        if ind == 1:
            cum_batch = x
        else:
            cum_batch = append_replicates(cum_batch, x)

        cum_batch.N_batches = n_batches
        cum_batch.gas_product = product
        acf_data = cum_batch.Compute_rate()
        print('Iteration ' + str(ind))
        print('Batches per trajectory: ' + str(cum_batch.Nbpt))
        print('Batch length (s): ' + str(cum_batch.batch_length))
        print('Rate: ' + str(cum_batch.rate))
        print('Rate confidence interval: ' + str(cum_batch.rate_CI))
        print('Autocorrelation: ' + str(cum_batch.ACF))
        print('Autocorrelation confidence: ' + str(cum_batch.ACF_CI))

        print('\n')

        batch_lengths.append(cum_batch.batch_length)
        rates_vec.append(cum_batch.rate)
        rates_vec_ci.append(cum_batch.rate_CI)
        acf_vec.append(cum_batch.ACF)
        acf_vec_ci.append(cum_batch.ACF_CI)

    print(batch_lengths)
    print(rates_vec)
    print(rates_vec_ci)
    print(acf_vec)
    print(acf_vec_ci)


    return cum_batch
