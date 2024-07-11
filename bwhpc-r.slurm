#!/bin/sh 
########## Begin MOAB/Slurm header ##########
#
#Partition
#SBATCH --partition=single

# Request number of nodes and CPU cores per node for job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# Estimated wallclock time for job
#SBATCH --time=00:10:00

########### End MOAB/SLURM header ##########

echo "Working Directory:                    $PWD"
echo "Running on host                       $HOSTNAME"
echo "Job id:                               $SLURM_JOB_ID"
echo "Job name:                             $SLURM_JOB_NAME"
echo "Number of nodes allocated to job:     $SLURM_JOB_NUM_NODES"
echo "Number of cores allocated to job:     $SLURM_JOB_CPUS_PER_NODE"


# Load module 
module load math/R/4.2.1

# Go to submit directory
cd $SLURM_SUBMIT_DIR
echo Directory is `pwd`

# Start R program
R CMD BATCH --no-save --no-restore fit.R

# notice for tests
echo "job terminated successfully"

