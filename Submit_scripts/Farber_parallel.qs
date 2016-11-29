#!/bin/bash
#$ -cwd
#$ -N parallel_python
#$ -pe mpi 40
#$ -l exclusive=1
#$ -o pp.out
#$ -e pp.err
#$ -l m_mem_free=2G
#$ -l h_cpu=24:00:00
#$ -m ea
# -M mpnunez@udel.edu

source /etc/profile.d/valet.sh
vpkg_rollback all
vpkg_require mpi,vtst

vpkg_require anaconda/2.5.0:python2
vpkg_require intel/2016
vpkg_require python-mpi4py/python2.7.8
vpkg_require python-numpy/python2.7.8
vpkg_require python-scipy/python2.7.8
vpkg_require python-matplotlib/python2.7.8

mpiexec -n 40 python FindSteadyState_AtoB.py