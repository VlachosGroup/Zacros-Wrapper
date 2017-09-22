# Read Zacros input files from a folder and then perform rate constant rescaling in a different folder

import zacros_wrapper as zw

''' ------------ User input section ------------ '''
exe_file = '/home/vlachos/mpnunez/bin/zacros_ZW.x'
KMC_source = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/stiff_input'
RunPath = '/home/vlachos/mpnunez/ZacrosWrapper/sample_systems/AtoB/ScaledownV3'
product_spec = 'B'                                  # product species
number_of_runs = 96
''' -------------------------------------------- '''

if __name__ == '__main__':                 # Make commands safe for MPI parallelization
	
    x = zw.kmc_traj(path = KMC_source)
    x.gas_prod = product_spec
    x.exe_file = exe_file
    x.ReadAllInput()
    final_data = zw.ReachSteadyStateAndRescale(x, RunPath, n_runs = 96, 
        n_batches = 1000, n_samples = 100, max_iterations = 10, parallel_mode = 'MPI')
