Python wrapper for the Zacros kinetic Monte Carlo (KMC) code
============================================================

This repository contains a Python library for a wrapper for the Zacros
kinetic Monte Carlo (KMC) code, which can be found at http://zacros.org/.
For users familiar with the Zacros software, our package offers ease-of-use
as well as additional analysis functionality. Modified source files for
Zacros are included which produce additional output files that are used by the wrapper.

* Documentation available at `<http://vlachosgroup.github.io/Zacros-Wrapper/>`_
* Download or clone source code from the  `Github repository <https://github.com/VlachosGroup/Zacros-Wrapper/>`_

Key features
------------
* Run KMC simulations with parallel processing
* Rescale rate constants of fast, equilibrated reactions to accelerate simulation
* Perform parameteric sensitivity analysis using the likelihood ratio method

Developers
----------
* Marcel Nunez (mpnunez28@gmail.com)
* Yifan Wang (wangyf@udel.edu)
* Taylor Robie
* Gerhard Wittreich, P.E.

Related Publications
---------------------
* `M. Nunez, T.A. Robie, D.G. Vlachos, J. Chem. Phys. 147, 164103 (2017). <http://aip.scitation.org/doi/full/10.1063/1.4998926>`_


Dependencies
-------------
* `Atomic simualtion environment <https://wiki.fysik.dtu.dk/ase/>`_ : Used to convert ab initio data to input parameters.
* `pMuTT <https://github.com/VlachosGroup/pMuTT/>`_ : Provides thermochemistry and kinetic parameters
* (optinoal) `mpi4py <http://pythonhosted.org/mpi4py/>`_ : Used for MPI parallelization of multiple trajectories. User can also use trivial parallelism with gridengine submit script.

Getting Started
----------------
1. Obtain modified Zacros executable. See Separate page.
2. Add the Zacros-Wrapper repository to your PYTHONPATH environment variable.
3. Configure input in demo file and run.
