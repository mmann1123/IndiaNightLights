#!/bin/bash

# set output and error output filenames, %j will be replaced by Slurm with the jobid
#SBATCH -o testing%j.out
#SBATCH -e testing%j.err 
 
# single node in the "short" partition
#SBATCH -N 1
#SBATCH -p short

# half hour timelimit
#SBATCH -t 2-00:00:00


#SBATCH --mail-type=ALL 
#SBATCH --mail-user=mmann1123@gwu.edu  

module load proj.4/4.8.0
module load gdal
module load R/3.0.2
module load gcc/4.9.0

# Run R code 
srun R CMD BATCH  ./grid_viirs_data_server_2015.R

