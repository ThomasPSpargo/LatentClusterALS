#!/bin/bash
#SBATCH --job-name=%paths%
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-%time%:00
#SBATCH --ntasks=8
#SBATCH --output=%logpath%%x.log

#Sort date into the variable $date, with the format YYYY_MM_DD
printf -v date '%(%Y_%m_%d)T' -1

#Create dated log directory if missing
mkdir -p %logpath%$date

#Run the R script for iterating Mplus model fits across files in the directory
Rscript /scratch/users/k1802739/ALS_LCA/scripts/01_Fitting_LCA/iterateModels.R \
	%dir%

#Once finished, generate fit statistics comparison
Rscript /scratch/users/k1802739/ALS_LCA/scripts/01_Fitting_LCA/graph_IC.R \
	%dir%bestout/

#Move log file to dated directory
mv %logpath%$SLURM_JOB_NAME.log \
	%logpath%$date/$SLURM_JOB_NAME.log