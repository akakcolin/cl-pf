#!/bin/bash
#
# modifed from humeniuka
# Run a command string on Slurm queue, for example
#
#   run_vasp.sh  "python run_vasp.py -i 1.vasp"  8 output
#

if [ $# == 0 ]
then
    echo " "
    echo "  Usage: $(basename $0)  command  nproc "
    echo " "
    echo "    submits a command (enclosed in quotation marks) to the queue with 'nproc' processors"
    echo " "
    echo "    in the current working directoy."
    echo " "
    echo "  Example:  $(basename $0)  'run_vasp.py -i POSCAR'  32  test"
    echo " "
    echo " "
    exit
fi

cmd=$1

# number of processors (defaults to 16)
nproc=${2:-16}

name=${3:-dft_task}
# errors and output of submit script will be written to this file
out=$(pwd)/${name}.out

# All options (arguments starting with --) are extracted from the command
# line and are passed on to sbatch.
options=""
for var in "$@"
do
    if [ "$(echo $var | grep "^--")" != "" ]
    then
  options="$options $var"
    fi
done

# The submit script is sent directly to stdin of qsub. Note
# that all '$' signs have to be escaped ('\$') inside the HERE-document.

# submit to PBS queue
#qsub <<EOF
# submit to slurm queue
sbatch $options <<EOF
#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --job-name=${name}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${nproc}

#SBATCH --output=${out}

module load compiler/2022.0.2
module load mkl/2022.0.2
module load mpi/2021.5.1
#module load cuda/11.8
export OMP_NUM_THREADS=1
# use ase
export VASP_PP_PATH=/data/psudopotential
export ASE_VASP_COMMAND="mpirun -np ${nproc} /data/bin/vasp_std"
export ASE_VASP_VDW=/data/bin/

#export LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH

# If the executable is in the submission directory, we want to be able
# to call it without the ./ prefix.
export PATH="\$SLURM_SUBMIT_DIR:\$PATH"

echo ------------------------------------------------------
echo SLURM_SUBMIT_HOST: \$SLURM_SUBMIT_HOST
echo SLURM_JOB_NAME: \$SLURM_JOB_NAME
echo SLURM_JOB_ID: \$SLURM_JOB_ID
echo SLURM_SUBMIT_DIR: \$SLURM_SUBMIT_DIR
echo SLURM_CPUS_ON_NODE: \$SLURM_CPUS_ON_NODE
echo ------------------------------------------------------
echo "Job is running on node(s):"
echo " \$SLURM_NODELIST "
echo ------------------------------------------------------
echo Path        : \$PATH
echo ------------------------------------------------------
echo Start date  : \$DATE
echo ------------------------------------------------------

# Here required modules are loaded and environment variables are set
export PYTHONUNBUFFERED=1

cd \$SLURM_SUBMIT_DIR

echo "Running command '$cmd'"
echo "   ##### START #####"
$cmd
echo "   ##### FINISH ####"

DATE=\$(date)
echo ------------------------------------------------------
echo End date: \$DATE
echo ------------------------------------------------------

EOF

echo "Submitting command '$cmd' (using $nproc processors)"
echo "Output will be written to '$out'."
