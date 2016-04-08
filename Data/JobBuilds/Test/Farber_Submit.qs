#!/bin/bash
#$ -cwd
#$ -pe mpi 1
#$ -l exclusive=0
#$ -l m_mem_free=2G
#$ -m as

source /etc/profile.d/valet.sh
#vpkg_rollback all
#export PYTHONPATH=/home/work/ccei_biomass/programs/ase/:$PYTHONPATH
#export PATH=/home/work/ccei_biomass/bin/:/home/work/ccei_biomass/programs/ase/tools:$PATH

#vpkg_require openmpi/1.8.2-intel64 python-numpy python-scipy

#export VASP_COMMAND="mpiexec -n 20 /home/work/ccei_biomass/programs/vasp5.4/vasp.5.4.1/build/std/vasp"
#export VASP_PP_PATH=/home/work/ccei_biomass/programs/vasp_psp/v54/


/home/1486/Zacros/ZacrosBinaryOutput_1.02/build/zacros.x
