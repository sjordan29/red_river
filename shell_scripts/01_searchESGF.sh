#!/bin/bash
#SBATCH -D /scratch/smj5vup/CMIP6/omoCMIP # working directory
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -t 0:30:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem-per-cpu=12288		# Memory per cpu (bytes)
#SBATCH --array=1-3			# Array of jobs to loop through
#SBATCH --mail-user=smj5vup@virginia.edu                                                                                                                           
#SBATCH --mail-type=ALL              

module load gcc
module load anaconda/2020.11-py3.8

# python scripts/data_setup/search_esgf_merge.py ${SLURM_ARRAY_TASK_ID}

python scripts/data_setup/test_separate_ds_search_esgf.py ${SLURM_ARRAY_TASK_ID}