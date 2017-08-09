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
	
    z = zw.RateRescaling()
    z.scale_parent_fldr = RunPath
    z.ReachSteadyStateAndRescale(product_spec, KMC_source, exe_file, n_runs = number_of_runs, platform = 'Squidward')