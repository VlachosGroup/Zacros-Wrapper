Examples
=========

We look at the examples in the 'examples' folder of the Github repository. For all scripts, the user must configure
the parameters (mostly directories) at the top of the script under the "User input" section.

Converting DFT data to Zacros input files
------------------------------------------------------

This script reads VASP data and computes parameters for the input files.

.. literalinclude:: ../examples/DFTtoThermoMain.py
   :language: python


Analyzing a single trajectory
-------------------------------

Single_trajectory.py analyzes a single Zacros trajectory and plots some data. 

.. literalinclude:: ../examples/Single_trajectory.py
   :language: python
    
This produces several .png images with data within the folder with the Zacros files.

Running and analyzing multiple trajectories
---------------------------------------------

Multiple_trajectories_MPI.py analyzes multiple trajectories, averages them, computes the rate, and performs sensitivity analysis. MPI is used to parallelize
running the trajectories. Therefore it must be launced with the ``mpiexec`` command.

.. literalinclude:: ../examples/Multiple_trajectories_MPI.py
   :language: python
    
This produces several .png files with statistics from the average of the trajectories.

Multiple_trajectories.py is usable on squidward.che.udel.edu and uses multiple job submissions to parallelize the jobs, as MPI does not allow for 
jobs spread over multiple nodes (i.e. > 16 cores).

.. literalinclude:: ../examples/Multiple_trajectories.py
   :language: python

Rate constant rescaling and steady state detection
-----------------------------------------------------

Scale_SS_Main.py takes Zacros input files and runs the statistical procedure described in our publication.
Parallelization is used to run multiple trajectories simultaneously. By default, it will use MPI parallelization, but 
the flag ``parallel_mode`` can be set to use the gridengine submission scripts on Squidward or Farber.

.. literalinclude:: ../examples/Rescale_and_steady_state.py
   :language: python