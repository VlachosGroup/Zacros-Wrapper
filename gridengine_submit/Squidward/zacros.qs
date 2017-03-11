#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe mpich 1
#


# The  executable:
export KMC_EXE="/home/vlachos/mpnunez/bin/zacros.x"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

${KMC_EXE}