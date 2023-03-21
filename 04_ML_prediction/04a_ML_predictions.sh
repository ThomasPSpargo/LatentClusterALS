#!/bin/bash

## This script prepares data for, and submits jobs to run, machine learning algorithms for predicting class membership.

######
### Setup
######

#Define the path containing the accepted Mplus model and appended cols
LCAfit=/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/savedata_C5/completesample/saved_mplusfit.tsv

#Define the path containing all scripts
scriptdir=/scratch/users/k1802739/ALS_LCA/scripts/04_ML_prediction

#Add a unique suffix for jobs submitted
suffix=_final

######
### Run
######
# Notes:
# This bash script first prepares the ML datasets from the results file passed to the DefineMLdatasets.R script.
# Then, across two loops, submit slutm jobs to run random forest (RF) and extreme gradient boosting (XGB) algorithms aiming to clasisfy each ML dataset.

#Prepare datasets 
Rscript ${scriptdir}/DefineMLdatasets.R --resultsfile $LCAfit
#Dataset files will be written to:
#ls $(dirname $LCAfit)/MLprediction/*_dataset.Rds

#Create a directory for the ML file logs
logdir=$(dirname $LCAfit)/MLprediction/logs
mkdir -p $logdir

printf "\n-------\n"

for file in $(dirname ${LCAfit})/MLprediction/*_dataset.Rds; do
	
	echo "Submit RF algorithm: $(basename $file)"
	
	# Submit Rscripts as batch jobs using a file generated within the printf call. This creates a slurm sbatch job file WITHOUT header and then declares settings within the sbatch call
	printf "#!/bin/bash\nRscript ${scriptdir}/random_forest_job.R --rdsDataset $file --tuneFor AUC --assignSuffix ${suffix} --ncores 10 --scriptDir ${scriptdir}" > ${scriptdir}/RF_jobFile.sh
		
	jobID=$(echo "RF_$(basename $file | sed -e s/_dataset.Rds//g)")
		
	sbatch --partition cpu --job-name $jobID --output ${logdir}/${jobID}.log --mem-per-cpu 10G --time 0-30:00 --ntasks 10 ${scriptdir}/RF_jobFile.sh
		
	#For organisational purposes, remove the file after submission
	rm ${scriptdir}/RF_jobFile.sh

done

printf "\n-------\n"

#Set the max number of cores to allow XGboost to use; set to 1 since script hangs when using parallel processing...
XGBncores=1

#Second loop across ML prepared datasets, submitting a series of XGB models partitioned by grid search
for file in $(dirname ${LCAfit})/MLprediction/*_dataset.Rds; do
	
	echo "Submit XGB: $(basename $file)"
	
	# Submit Rscripts as batch jobs using a file generated within the printf call. This creates a slurm sbatch job file WITHOUT header and then declares settings within the sbatch call
	printf "#!/bin/bash\nRscript ${scriptdir}/xgboost_stagedTuning.R --rdsDataset ${file} --tuneFor AUC --assignSuffix ${suffix} --scriptDir ${scriptdir} --ncores $XGBncores --tuningLogFile tuneinfo.log" > ${scriptdir}/XGB_staged_jobFile.sh
		
	jobID=$(echo "XGB_$(basename $file | sed -e s/_dataset.Rds//g)${suffix}")
		
	sbatch --partition cpu --job-name $jobID --output ${logdir}/${jobID}.log --mem-per-cpu 10G --time 0-40:00 --ntasks $XGBncores ${scriptdir}/XGB_staged_jobFile.sh
		
	#For organisational purposes, remove the file after submission
	rm ${scriptdir}/XGB_staged_jobFile.sh

done
		