#!/bin/bash
#SBATCH -D /scratch/smj5vup/CMIP6/omoCMIP/ # working directory
#SBATCH -o /scratch/smj5vup/CMIP6/output/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH -N 2				# number of nodes 
#SBATCH --ntasks-per-node 40
#SBATCH -p parallel        				# Queue name "dev"     					# allocation name
#SBATCH -t 24:00:00       					# Run time (hh:mm:ss) - up to 36 hours // 1 hr for dev queue
#SBATCH --mail-user=smj5vup@virginia.edu      # address for email notification
#SBATCH --mail-type=ALL  


module load gcc
module load openmpi
module load anaconda/2020.11-py3.8

srun python scripts/downscale/mpi_downscale_tair.py

# python scripts/downscale/SJ_test_pr.py 1