#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

#$ -N zacros_JA 					#This is the name of the job array
#$ -t 1-4  							#Assumes task IDs increment by 1; can also increment by another value
#$ -tc 4 							#This is the total number of tasks to run at any given moment
#$ -pe openmpi-smp 1 				#Change the last field to the number of processors desired per task

job_file='./dir_list.txt'
#Change to the job directory
job_path=$(sed -n "$SGE_TASK_ID p" "$job_file")
cd "$job_path" #SGE_TASK_ID is the task number in the range <task_start_index> to <task_stop_index>
                  #This could easily be modified to take a prefix; ask me how.

# The  executable:
export KMC_EXE="/home/vlachos/mpnunez/bin/zacros.x"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

${KMC_EXE}