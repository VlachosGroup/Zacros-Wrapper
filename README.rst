Python wrapper for the Zacros kinetic Monte Carlo (KMC) code
============================================================

This repository contains a Python library for a wrapper for the Zacros 
kinetic Monte Carlo (KMC) code, which can be found at http://zacros.org/. 
For users familiar with the Zacros software, our package offers ease-of-use 
as well as additional analysis functionality. Modified source files for 
Zacros are included which produce additional output files that are used by the wrapper.

* Documentation available at `<http://vlachosgroup.github.io/Zacros-Wrapper/>`_ (under development)
* Download or clone source code from the  `Github repository <https://github.com/VlachosGroup/Zacros-Wrapper/>`_

Key features
------------
* Run KMC simulations with parallel processing (uses `mpi4py <http://pythonhosted.org/mpi4py/>`_ or job submission scripts for University of Delaware clusters)
* Rescale rate constants of fast, equilibrated reactions to accelerate simulation
* Perform parameteric sensitivity analysis using the likelihood ratio method

Developers
----------
* Marcel Nunez (mpnunez@udel.edu)
* Taylor Robie
* Gerhard Wittreich, P.E.

Related publications
-----------------------
* M. Nunez, T.A. Robie, and D.G Vlachos, “Acceleration and Sensitivity Analysis of Lattice Kinetic Monte Carlo Simulations Using Parallel Processing and Rate Constant Rescaling” (in preparation)

Getting Started
----------------
1. Obtain modified Zacros executable. See Separate page.
2. Add the Zacros-Wrapper repository to your Python path.
3. Configure input in demo file and run.