#!/bin/bash
#$ -cwd
#$ -j y
# -N ZW_COox_FDSA
#$ -S /bin/bash
#$ -pe openmpi-smp 1
#

# Get our environment setup:
vpkg_require "python/2.7.8"
vpkg_require "openmpi/1.6.3-gcc"
vpkg_require "python-numpy"
vpkg_require "python-scipy"
vpkg_require "python-mpi4py"

# The  executable:
export PYTHON_EXE="python Run_replicates.py"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

${PYTHON_EXE}