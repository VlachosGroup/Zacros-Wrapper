#!/bin/bash
#$ -cwd
#$ -j y
#$ -N mpi_test
#$ -S /bin/bash
#$ -pe openmpi-smp 16
#

# Get our environment setup:
vpkg_require "python/2.7.8"
vpkg_require "openmpi/1.6.3-gcc"
vpkg_require "python-numpy"
vpkg_require "python-scipy"
vpkg_require "python-mpi4py"

# The  executable:
mpiexec -n 16 python FindSteadyState.py