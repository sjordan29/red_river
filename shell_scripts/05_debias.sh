#!/bin/bash
#SBATCH -D /scratch/smj5vup/CMIP6/omoCMIP/ # working directory
#SBATCH -o /scratch/smj5vup/CMIP6/output/job.%j.%N.out   # Name of the output file (eg. myMPI.oJobID)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 20
#SBATCH -p standard          				# Queue name "dev"
#SBATCH -t 24:00:00       					# Run time (hh:mm:ss) - up to 36 hours // 1 hr for dev queue 
#SBATCH --mail-user=smj5vup@virginia.edu      # address for email notification
#SBATCH --mail-type=ALL                  	# email at Begin and End of job

module purge 

# r dependencies 
module load gcc/7.1.0  
module load openmpi/3.1.4
module load intel/20.0  
module load intelmpi/20.0
module load goolf/7.1.0_3.1.4 R/3.6.3 
module load HDF5/1.10.4

# package dependencies 
module load hdf5/1.10.6 
module load netcdf/4.7.3

Rscript scripts/bias_correction/debias_gcm_data.R ${SLURM_CPUS_PER_TASK}