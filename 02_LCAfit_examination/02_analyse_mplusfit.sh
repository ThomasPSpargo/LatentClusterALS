#!/bin/bash
#SBATCH --job-name=analyse_mplus_model
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0-1:00
#SBATCH --ntasks=1
#SBATCH --output=/scratch/users/k1802739/ALS_LCA/logs/%x.log

## This script submits a series of (mostly) R scripts to perform various analyses, comparing classes determined by Mplus

######
### Setup
######

#Sort date into the variable $date, with the format YYYY_MM_DD
printf -v date '%(%Y_%m_%d)T' -1

#Use this directory for all log files
logpath=/scratch/users/k1802739/ALS_LCA/logs

#Create dated log directory if missing
mkdir -p $logpath/$date

#Set path to directory containing all scripts
scriptpath="/scratch/users/k1802739/ALS_LCA/scripts/02_LCAfit_examination"

# Set path to dataset containing all data fields to be matched with the subset used by  MPLUS for LCA
rootdataset="/scratch/users/k1802739/ALS_LCA/final_datasets/JNT_LCA_Allfields.tsv"

# Specify the column storing numeric IDs that were passed to mplus (note that this column will be called ID in the Mplus output scripts themselves)
rootIDcolumn="new_num_id"

#Predefine this option
#--MplusKeepcolumns a regex pattern passed to grepl which declares columns from the Mplus output file that should be retained when combining datasets.  the `|` indicates 'or' and this string matches columns "C", "CPROB1"..."CPROBk" where k is the number of classes.
mplusKeepcolumns="^C\$|^CPROB"

#Define column names that will define the clinical parameters used within the LCA model
#Some of these columns are generated during the 'process_dataset.R' script
predcols="CTYDEL,AGEONS_nml,SRV_nml,Site_of_Onset,Sex_at_birth,Phenotype"

#Path to the secondary dataset (used for applying PM discovery data to STRENGTH validation data)
validationDataDir="/scratch/users/k1802739/ALS_LCA/final_datasets/mplusfiles"

#Path to mplus helper functions directory
helperFunDir=/scratch/users/k1802739/ALS_LCA/scripts/mplusHelper


#Specify the file paths to the LCA 5-class models for the Joint (PM+STRENGTH) dataset, and for PM only. Analyse both the full data and the data which omits people with missing disease duration or diagnostic delay
declare -a PATHS=("/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/JNT_FULL_ANY_c5_re3200_stscale_15_trydbl.out" "/scratch/users/k1802739/ALS_LCA/FinalModels/PM_FULL_ANY/bestout/PM_FULL_ANY_c5_re3200_stscale_15_trydbl.out" "/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_DDOMIT/bestout/JNT_FULL_DDOMIT_c5_re3200_stscale_15_trydbl.out" "/scratch/users/k1802739/ALS_LCA/FinalModels/PM_FULL_DDOMIT/bestout/PM_FULL_DDOMIT_c5_re3200_stscale_15_trydbl.out")

######
### Begin
######

#Loop across files in the array
for FILE in ${PATHS[@]}; do
	
	printf "##################################\n---------BEGIN NEXT LOOP----------\n##################################\n"
	echo "This is the best fit model $FILE"
	
	printf "\nExtracting the model reappending non-model columns\n"

	### Generate model data save files and store the output to the file 'saved_mplusfit.tsv'
	#--Mplusoutfile identifies the .out file specified
	#--PreMplusFields identifies the dataset from which the records used in mplus were derived
	#--PreMplusIDname names the ID column against which the Mplus outputs should be matched
	#--MplusKeepcolumns regex pattern declaring columns from the Mplus output file that should be retained when combining datasets.  the `|` indicates 'or' and this string matches columns "C", "CPROB1"..."CPROBK" where k is the number of classes.
	Rscript ${scriptpath}/obtainSavedata.R \
		--Mplusoutfile $FILE \
		--PreMplusFields $rootdataset \
		--PreMplusIDname  $rootIDcolumn \
		--MplusKeepcolumns $mplusKeepcolumns \
		--helperFunDir $helperFunDir

	#Determine the number of classes in the accepted model
	nclasses=$(basename $FILE | gawk 'match($0, /.*(c[0-9+]).*/,a) {print a[1]}' | tr '[:lower:]' '[:upper:]')
	
	#Declare the output file from the obtainSavedata script
	initresultfile="$(dirname $FILE)/savedata_${nclasses}/saved_mplusfit.tsv"
	
	### Rewrite the file including normalisation columns.
	# Subset additionally into a complete dataset and a 'censoring-excluded' cohort
	
	#Mplus class order is arbitrary, if this is TRUE class numbers will be reordered by sample size largest to smallest.
	#Note that the SVALUES argument of mplus should handle class reordering. However, mplus support confirmed our report of a software bug that means the functionality does not work properly with the time-to-event/censored variable; people's classifications were changed when they should not be.
	reorderClasses=TRUE
	
	printf "\n---------------------\nParse the output and prepare for analysis, recode the class order? ${reorderClasses}\n"

	#--filetoparse declares the input file
	#--writePanel TRUE #Stores a standardisation panel for continuous variables, so that equivalent normalisation can be applied to secondary data
	#--writeSrvSubset TRUE #Setting TRUE dictates that files should be returned in a subdirectory, splitting between the complete sample and a subsample of people with non-censored survival
	Rscript ${scriptpath}/process_dataset.R \
		--filetoparse $initresultfile \
		--writePanel TRUE \
		--writeSrvSubset TRUE \
		--reorderClasses $reorderClasses
		
	#Run all steps for both the complete and noncensored versions of the dataset
	declare -a resultversions=("completesample" "noncensored")
	for version in ${resultversions[@]}; do
		
		#Declare the results file to use
		resultfile="$(dirname $FILE)/savedata_${nclasses}/${version}/saved_mplusfit.tsv"

		#Apply r script summary functions to the final csv file
		printf "\n##################\n Running scripts for the file: $resultfile"
		
		printf "\n---------------------\nGenerate descriptive statistics\n"

		Rscript ${scriptpath}/generate_descriptives.R --resultsfile $resultfile

		printf "\n---------------------\nSummarise class probabilities\n"

		#Pass the dataset to summarise average class probabilities in table and figure
		#This looks towards the cprobs columns
		Rscript ${scriptpath}/ave_cprobs.R $resultfile

		printf "\n---------------------\nApply survival analyses\n"
		
		# Run survival analysis - note that predictor columns are defined WITHIN the script here so that their names can be updated
		#--keepcols specify the columns to retain in the dataset
		#--coeftableFunpath points directly to custom function for preparing a pretty coefficients table and forest plotting
		Rscript ${scriptpath}/survival_analysis.R \
			--resultsfile $resultfile \
			--coeftableFunpath "${scriptpath}/coxCoefTable.R"
		
		printf "\n---------------------\nRun linear discriminant analysis\n"

		#--keepcols specify the columns to retain in the dataset
		#--Predictorcols indicates the columns to be used as predictors
		#--intPredict Runs LDA with an internal prediction step, after running on the full sample, running again after splitting into 80:20 train and test
		Rscript ${scriptpath}/linear_discriminant_analysis.R \
			--resultsfile $resultfile \
			--keepcols "C,CTYDEL,AGEONS_nml,SRV_nml,Site_of_Onset,Sex_at_birth,Phenotype,origin" \
			--Predictorcols $predcols \
			--intPredict TRUE
		
		
		printf "\n---------------------\nRun multinomial logistic regression\n"

		Rscript ${scriptpath}/Mnom_Logistic_regression.R \
			--resultsfile $resultfile \
			--keepcols "C,CTYDEL,AGEONS_nml,SRV_nml,Site_of_Onset,Sex_at_birth,Phenotype,origin" \
			--Predictorcols $predcols

		#When analysing the discovery sample/full data, apply the fit to secondary datasets
		if [[ $FILE =~ "PM" &&  $version =~ "completesample" ]]; then
			printf "\n---------------------\nApply the model to new data\n"

			#Identify the STRENGTH dataset to compare against (i.e ANY or DDOMIT)
			ptype=$(basename ${FILE} | cut -d'.' -f 1 | gawk 'match($0, /.*_([A-z]+)_c[0-9]+.*/,a) {print a[1]}')
			DSET=$validationDataDir/STR_NOPM_${ptype}.csv
				
			#Construct path to and prefix for output filepath for secondary data fits
			directory=$(dirname $resultfile)
			prefix=$(basename $DSET | sed 's/\..*//g')
			
			#Identify the mplus output file describing the model fit
			savefile=$(echo $(dirname $initresultfile)/*.out | xargs ls)
			
			printf "##################################\nPREDICT IN $DSET\n##################################\n"

			#--DEFINE and --USEOBSERVATIONS options can be used to pass any filtering instructions to MPLUS

			#Repeat modelling in secondary dataset 
			predictfile=$(Rscript ${scriptpath}/predictNewData.R \
							--SOURCEFILE $savefile \
							--DATASET $DSET \
							--out "${directory}/${prefix}" \
							--PreMplusFields $rootdataset \
							--PreMplusIDname $rootIDcolumn \
							--MplusKeepcolumns $mplusKeepcolumns \
							--helperFunDir $helperFunDir)	


			printf "\n---------------------\nParse the predicted file, standardising the continuous variables relative to discovery dataset\n"

			#--normalisationPanel indicates where the normalisation panel for continuous variables is saved
			#See comments above describing previous application	
			Rscript ${scriptpath}/process_dataset.R \
			--filetoparse $predictfile \
			--normalisationPanel "$(dirname $initresultfile)/standardisationPanel.tsv" \
			--reorderClasses $reorderClasses

			printf "\n---------------------\nGenerate descriptive statistics for secondary dataset\n"

			Rscript ${scriptpath}/generate_descriptives.R --resultsfile $predictfile

			
			printf "\n---------------------\nSummarise class probabilities for new data\n"

			#Generate average cprobabilities plot
			Rscript ${scriptpath}/ave_cprobs.R $predictfile

			printf "\n---------------------\nApply K-Nearest Neighbours to test independent classification approach\n"

			#--resultsfile indicates the DISCOVERY dataset
			#--dataset2 indicates the validation dataset, used as 'test' data in KNN
			Rscript ${scriptpath}/KNNvalidation.R \
				--resultsfile $resultfile \
				--keepcols "C,CTYDEL,AGEONS_nml,SRV_nml,Site_of_Onset,Sex_at_birth,Phenotype,origin" \
				--Predictorcols $predcols \
				--dataset2 $predictfile	
			
		fi
	done
done

#Move log file to dated directory
mv $logpath/$SLURM_JOB_NAME.log \
	$logpath/$date/$SLURM_JOB_NAME.log