#!/bin/bash
#SBATCH -D /scratch/smj5vup/CMIP6/omoCMIP/	 # working directory
#SBATCH -o /scratch/smj5vup/CMIP6/output/job.%j.%N.out	# Name of the output file (eg. myMPI.oJobID)
#SBATCH --ntasks=1				# Number of tasks per serial job (must be 1)
#SBATCH -p standard				# Queue name "standard" (serial)
#SBATCH -t 5:00:00				# Run time per serial job (hh:mm:ss)
#SBATCH --mem-per-cpu=12288		# Memory per cpu (bytes)
#SBATCH --array=2-50			# Array of jobs to loop through
#SBATCH --mail-user=smj5vup@virginia.edu                                                                                                                           
#SBATCH --mail-type=ALL    

module load gcc/7.1.0  
module load openmpi/3.1.4
module load intel/20.0  
module load intelmpi/20.0
module load goolf/7.1.0_3.1.4 R/3.6.3 
module load gcc/7.1.0 mvapich2/2.3.1
module load gcc/7.1.0  openmpi/3.1.4
module load intel/18.0  mvapich2/2.3.3
module load grass

# hdf dependencies 
module load hdf5/1.10.6


Rscript scripts/swat_setup/precip_zonal_stat.R ${SLURM_ARRAY_TASK_ID}