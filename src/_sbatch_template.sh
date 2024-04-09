#!/bin/bash

#SBATCH -J _jobname 
#SBATCH -o _jobnamelogfile
#SBATCH -e _jobnameerrfile
#SBATCH -N 1
#SBATCH -n 36

#SBATCH -p  test


source /data/intel/parallel_studio_xe_2019/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_FC=ifort
export I_MPI_F90=ifort
export PATH=/bin/:$PATH
export OMP_NUM_THREADS=1
ulimit -s unlimited
#_file
_input
