Installation instructions
-------------------------
* Download the Zacros files from the [Zacros website](http://www.e-lucid.com/i/software/Zacros.html).
* Replace the .f90 files in main folder with the ones in the Zacros_mod folder from the ZacrosWrapper repository.
* Compile with cmake following the compiling instructions provided with Zacros. Using cmake ensures that the executable will run on the platform on which it is compiled.
* Rename the executable file (.x in Unix, .exe in Windows) and put it in a directory that the Wrapper can access.

Additional output files
-----------------------

+------------------------+--------------------------------------------------+
| file name              | description                                      |
+========================+==================================================+
| PropCounter_output.bin | integral propensities, used for better averaging |
+------------------------+--------------------------------------------------+
| Prop_output.bin        | reaction propensities                            |
+------------------------+--------------------------------------------------+
| Hist.bin               | lattice states                                   |
+------------------------+--------------------------------------------------+
| E.bin                  | lattice energy                                   |
+------------------------+--------------------------------------------------+
| clusterocc.bin         | cluster counts                                   |
+------------------------+--------------------------------------------------+
| SA_output.bin          | trajectory derivatives                           |
+------------------------+--------------------------------------------------+
| specnum_output.bin     | surface species populations                      |
+------------------------+--------------------------------------------------+

Compiling on [Farber](http://farber.hpc.udel.edu/) (University of Delaware computing cluster)
---------------------------------------------------------------------------------------------

If you wish to compile Zacros on Farber you will need to make some configurational changes as the default compile settings are very inefficient.

* Import the intel compiler: ``vpkg_require intel/2016``
* ``export CMAKE_Fortran_COMPILER=ifort``
* Usual installation

``
cd path/to/source/of/Zacros 
mkdir build` # If it does not exist yet   
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
``
* Instead of ``make`` as usual, type:  ``ccmake`` .
* Once within the ccmake, press ``t`` for advanced mode, then set the following flags
``
CMAKE_BUILD_TYPE=Release
OMP_NUM_PROCS=1
OpenMP_Fortran_FLAGS=-qopenmp
CMAKE_Fortran_COMPILER=ifort
``
* Press ``c`` to configure, and then ``g`` to generate and exit. Finally: ``make``  