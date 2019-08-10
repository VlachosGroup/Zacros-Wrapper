# Read the output for a Zacros trajectory and plot data

import zacros_wrapper as zw
import os

''' ------------ User input section ------------ '''
RunPath = os.path.join(os.path.dirname(__file__), '../input_files/AtoB')

''' -------------------------------------------- '''

''' Read simulation results '''
my_trajectory = zw.kmc_traj()                           # Create single trajectory object
my_trajectory.Path = RunPath                            # Set directory for files
my_trajectory.ReadAllOutput(build_lattice=True)         # Read input and output files

''' Plot data '''
my_trajectory.PlotSurfSpecVsTime()                      # Plot surface species populations versus time
my_trajectory.PlotGasSpecVsTime()                       # Plot gas species populations versus time
my_trajectory.PlotElemStepFreqs(time_norm = True)       # Plot elementary step frequencies
my_trajectory.PlotLattice()                             # Plot the lattice
my_trajectory.LatticeMovie()                            # Plot the lattice snapshots
