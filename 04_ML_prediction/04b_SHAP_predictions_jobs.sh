#!/bin/bash

######
### Setup
######

#Declare outbound filepath and prefix ('shapVals') RDS files with shap values for each dataset/algorithm
outfile=/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/savedata_C5/completesample/MLprediction/weightedShapley/shapVals
mkdir -p $(dirname $outfile)

#Identify ML to extract from by supplying filepaths one per row in a text file
fits=/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/savedata_C5/completesample/MLprediction/MLfitsForShapley.txt

#Delare path for shap job script
scriptdir=/scratch/users/k1802739/ALS_LCA/scripts/04_ML_prediction

#Declare directory containing the datasets upon which models were trained
datasetdir=/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/savedata_C5/completesample/MLprediction

#Declare the directory in which to return logs for the submitted slurm jobs
logdir=/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/savedata_C5/completesample/MLprediction/logs

#nPredictSample=NA #Indicate this argument to predict in a data subsample

#For each ML model prepare a batch job file and submit 
loops=$(cat $fits | wc -l) 	

for i in $(seq 1 $loops); do	#Loop for each machine learning fit

	file=$(awk 'NR=='$i' {print}' $fits)
	#--nPredictSample ${nPredictSample} Exclude this argument for the full fit
	printf "#!/bin/bash\nRscript ${scriptdir}/SHAP_job.R --finalTune ${file} --datasetDir ${datasetdir} --ntasks 1 --supplyWeights TRUE --outfile ${outfile}" > ${scriptdir}/SHAP_jobFile.sh
		
	jobID="SHAP_file_${i}"
				
	sbatch --partition cpu --job-name $jobID --output ${logdir}/${jobID}.log --mem-per-cpu 10G --time 0-20:00 --ntasks 1 ${scriptdir}/SHAP_jobFile.sh
	
	#For organisational purposes, remove the file after submission
	rm ${scriptdir}/SHAP_jobFile.sh

done