#!/bin/bash
#$ -cwd
#$ -pe mpi 1
#$ -l exclusive=0
#$ -l m_mem_free=2G
#$ -m as

source /etc/profile.d/valet.sh
vpkg_require openmpi/1.10.2-intel64-2016

/home/1486/Zacros/ZacrosBinaryOutput_1.02/build/zacros.x
