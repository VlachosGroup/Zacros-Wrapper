#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

#$ -N zacros_JA 					#This is the name of the job array
#$ -t 1-20  							#Assumes task IDs increment by 1; can also increment by another value
#$ -tc 10 							#This is the total number of tasks to run at any given moment
#$ -pe threads 1 				#Change the last field to the number of processors desired per task
#
# Change the following to #$ and set the amount of memory you need
# per-slot if you're getting out-of-memory errors using the
# default:
#$ -l m_mem_free=4G

source /etc/profile.d/valet.sh

# Use vpkg_require to setup the environment:
vpkg_require intel/2016

# Ensure that the OpenMP runtime knows how many processors to use;
# Grid Engine automatically sets NSLOTS to the number of cores granted
# to this job:
export OMP_NUM_THREADS=$NSLOTS

job_file='./dir_list.txt'
#Change to the job directory
job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")
cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>
                  #This could easily be modified to take a prefix; ask me how.

# Now append whatever commands you use to run your OpenMP code:
/home/1483/bin/zacros_ZW.x