Examples
=========

We look at the examples in the 'examples' folder of the Github repository.

Converting DFT data to Zacros input files
------------------------------------------------------

This script analyzes a single Zacros trajectory and plots some data. 

.. literalinclude:: ../examples/DFTtoThermoMain.py
   :language: python
    
This produces ...


Analyzing a single trajectory (Single_trajectory.py)
------------------------------------------------------

This script analyzes a single Zacros trajectory and plots some data. 

.. literalinclude:: ../examples/Single_trajectory.py
   :language: python
    
This produces ...

Analyzing multiple trajectories (Multiple_trajectories.py)
-----------------------------------------------------------

This script analyzes multiple trajectories, averages them, computes the rate, and performs sensitivity analysis.

.. literalinclude:: ../examples/Multiple_trajectories.py
   :language: python
    
This produces ...

Rate constant rescaling and steady state detection (Scale_SS_Main.py)
-----------------------------------------------------------------------

This script takes Zacros input files and runs the statistical procedure described in our publication.
Parallelization is used to run multiple trajectories simultaneously.

.. literalinclude:: ../examples/Scale_SS_Main.py
   :language: python
    
This produces ...