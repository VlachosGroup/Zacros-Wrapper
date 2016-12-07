# -*- coding: utf-8 -*-
"""
Created on Mon Dec 05 12:40:40 2016

@author: mpnun
"""

# Read a simulation and set up the next iteration of a scaledown procedure
# INCOMPLETE - Trying to make job work on Farber without having to resort to manual tweaks

import os
import sys
import copy
import numpy as np

sys.path.append('/home/1483/ZacrosWrapper')
import ZacrosWrapperPy as zw

scale_parent_fldr = '/home/1483/ZacrosWrapper/KMC_data/WGS/Scaledown/'

include_stiff_reduc = True
iteration = 1
n_runs = 1000
Product = 'CO2'
ss_inc = 2.0

''' Read and analyze current iteration '''

cur_batch = zw.Replicates()
iter_fldr = scale_parent_fldr + 'Iteration_' + str(iteration) + '/'
cur_batch.ParentFolder = iter_fldr
cur_batch.ReadMultipleRuns()
SDF_vec = np.ones(cur_batch.runtemplate.Reactions['nrxns'])

# Add data to running list
#if iteration == 1:
#    pass            # Do not use data from first run because it is on event rather than time
#elif iteration == 2:
#    cum_batch = copy.deepcopy(cur_batch)
#elif iteration > 2:
#    for run_ind in range(n_runs):
#        cum_batch.runList[run_ind] = KMC_Run.time_sandwich(cum_batch.runList[run_ind], cur_batch.runList[run_ind])

cum_batch = copy.deepcopy(cur_batch)

# Test steady-state
if iteration == 1:
    is_steady_state = False
else:
    cum_batch.AverageRuns()
    cum_batch.ParentFolder = iter_fldr
    cum_batch.runAvg.Path = iter_fldr
    correl = cum_batch.CheckAutocorrelation(Product)
    not_change = cum_batch.runAvg.CheckSteadyState(Product)
    is_steady_state = np.abs(correl[0]) < 0.05 and not_change
    
    # Record information about the iteration
    cum_batch.runAvg.CalcRateTraj(Product)

    cum_batch.runAvg.PlotSurfSpecVsTime()        
    cum_batch.runAvg.PlotIntPropsVsTime()
    cum_batch.runAvg.PlotRateVsTime()  

# Test stiffness
cur_batch.AverageRuns()
scaledown_data = zw.Replicate_Analysis.ProcessStepFreqs(cur_batch.runAvg)         # compute change in scaledown factors based on simulation result
delta_sdf = scaledown_data['delta_sdf']
rxn_speeds = scaledown_data['rxn_speeds']
if include_stiff_reduc:
    stiff_cutoff = 1
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


''' Prepare and run next iteration '''


# Make folder for iteration
iter_fldr = scale_parent_fldr + 'Iteration_' + str(iteration) + '/'
if not os.path.exists(iter_fldr):
    os.makedirs(iter_fldr)
    
# Create object for batch
cur_batch = zw.Replicates()
cur_batch.ParentFolder = iter_fldr
cur_batch.n_runs = n_runs
cur_batch.Product = Product

cur_batch.runtemplate = zw.KMC_Run()
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
    cur_batch.runtemplate.Report['hist'] = ['event', max_events]       # only record the initial and final states

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