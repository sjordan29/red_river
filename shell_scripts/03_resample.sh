#!/bin/bash
#SBATCH -D /scratch/smj5vup/CMIP6/omoCMIP # working directory
#SBATCH -o /scratch/smj5vup/CMIP6/output/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -t 0:30:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem-per-cpu=12288		# Memory per cpu (bytes)
#SBATCH --array=1-1			# Array of jobs to loop through
#SBATCH --mail-user=smj5vup@virginia.edu                                                                                                                           
#SBATCH --mail-type=ALL              


# dependencies
module load gcc/9.2.0
module load netcdf/4.7.3
module load hdf5/1.10.6
module load intel/18.0
module load jasper/1.900.1
module load grib_api/1.28.0

# cdo - local 
module use $HOME/modulefiles
module load cdo 

# anaconda 
module load anaconda/2020.11-py3.8


python scripts/data_setup/04_resample_cmip6.py

