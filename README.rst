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
* Run KMC simulations with parallel processing (uses GridEngine job arrays)
* Rescale rate constants of fast, equilibrated reactions to accelerate simulation
* Perform parameteric sensitivity analysis using finite-difference or likelihood ratio methods

Developers
----------
* Marcel Nunez (mpnunez@udel.edu)
* Taylor Robie
* Gerhard Wittreich, P.E.

Related publications
-----------------------
* M. Nunez, T.A. Robie, and D.G Vlachos, “Acceleration and Sensitivity Analysis of Lattice Kinetic Monte Carlo Simulations Using Parallel Processing and Rate Constant Rescaling” (in preparation)