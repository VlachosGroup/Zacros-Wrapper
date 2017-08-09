# Read the output for a Zacros trajectory and plot data

import zacros_wrapper as zw

''' ------------ User input section ------------ '''
RunPath = '/home/vlachos/wangyf/Alumina/mechanismI/5'
''' -------------------------------------------- '''

''' Read simulation results '''
my_trajectory = zw.kmc_traj()                           # Create single trajectory object
my_trajectory.Path = RunPath                            # Set directory for files
my_trajectory.ReadAllOutput(build_lattice=True)         # Read input and output files

''' Plot data '''
my_trajectory.PlotSurfSpecVsTime()                                      # Plot surface species populations versus time
my_trajectory.PlotGasSpecVsTime()                                       # Plot gas species populations versus time
my_trajectory.PlotElemStepFreqs(site_norm = n_Pd, time_norm = True)     # Plot elementary step frequencies
my_trajectory.PlotLattice()                                             # Plot the lattice
my_trajectory.LatticeMovie()                                            # Plot the lattice snapshots