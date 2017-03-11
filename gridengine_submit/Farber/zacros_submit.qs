#
# Template:  Generic OpenMP Jobs
# Revision:  $Id: openmp.qs 523 2014-09-16 14:29:54Z frey $
#
# Usage:
# 1. Modify "NPROC" in the -pe line to reflect the number
#    of processors desired
# 2. Jobs default to using 1 GB of system memory per slot.  If you
#    need more than that, set the m_mem_free complex.
#
#$ -pe threads 1
#
# Change the following to #$ and set the amount of memory you need
# per-slot if you're getting out-of-memory errors using the
# default:
#$ -l m_mem_free=4G
#
# If you want an email message to be sent to you when your job ultimately
# finishes, edit the -M line to have your email address and change the
# next two lines to start with #$ instead of just #
# -m eas
# -M my_address@mail.server.com
#

source /etc/profile.d/valet.sh

# Use vpkg_require to setup the environment:
vpkg_require intel/2016

# Ensure that the OpenMP runtime knows how many processors to use;
# Grid Engine automatically sets NSLOTS to the number of cores granted
# to this job:
export OMP_NUM_THREADS=$NSLOTS

#
# Now append whatever commands you use to run your OpenMP code:
#
/home/1486/Zacros/ZacrosBinaryOutput_1.02/build-jeff/zacros.x

