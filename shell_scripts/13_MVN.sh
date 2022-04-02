#!/bin/bash
#SBATCH -D /scratch/smj5vup/CMIP6/omoCMIP # working directory
#SBATCH -o /scratch/smj5vup/CMIP6/output/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -t 3:30:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem-per-cpu=12288		# Memory per cpu (bytes)
#SBATCH --array=1-1			# Array of jobs to loop through
#SBATCH --mail-user=smj5vup@virginia.edu                                                                                                                           
#SBATCH --mail-type=ALL              


module load anaconda/2020.11-py3.8

python scripts/trend_analysis/06_MVN.py