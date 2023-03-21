#!/bin/bash
#SBATCH --job-name=Genetic_comparisons
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-1:00
#SBATCH --ntasks=1
#SBATCH --output=/scratch/users/k1802739/ALS_LCA/logs/%x.log


######
### Setup 	
######

#Specify the file from the accepted dataset
LCAfit=/scratch/users/k1802739/ALS_LCA/FinalModels/JNT_FULL_ANY/bestout/savedata_C5/completesample/saved_mplusfit.tsv

#Set path to directory containing all scripts
scriptpath=/scratch/users/k1802739/ALS_LCA/scripts/03_Genetic_analyses

#Declare the control dataset location
controlData=/scratch/users/k1802739/ALS_LCA/final_datasets/PM_CONTROLS_data.tsv

#Declare column names for rare genetic and PRS features in comma separated list
raregenetics="SOD1,c9_status,AutoProt,RNAFunc,CytoTransp,BurdenG"

prscols="ALScc_SBayesR,FTDcc_SBayesR,SZcc_SBayesR,AZcc_SBayesR,PDcc_SBayesR"

######
### RUN
######
# Notes:
# The script runs the fishers exact test script 1 time only
# The analysis with PRS are performed 4 times:
# 	1. a simple model comparing a given class to all other classes
# 	2. a model with a given PRS and principal components 1-5, comparing a given class to all other classes
# 	3. a simple model comparing a given class to the healthy control cohort
# 	4. a model with a given PRS and principal components 1-5, comparing a given class to the healthy controls
# A final script is run that visualises the results of all PRS analyses

	
start=`date +%s`
	
	printf "##########\nRun the Fisher's exact test script\n##########\n"	
	
	Rscript ${scriptpath}/fishers_compare.R \
		--resultsfile $LCAfit \
		--Predictorcols $raregenetics
		
	printf "##########\nRun the PRS script\n##########\n"	
		
	PRSsimple=$(Rscript	${scriptpath}/prs_compare.R \
		--resultsfile $LCAfit \
		--Predictorcols $prscols \
		--IncludeNofPCs 0)
		
	printf "##########\nRun the PRS script with covariates\n##########\n"	
		
	PRScovar=$(Rscript	${scriptpath}/prs_compare.R \
		--resultsfile $LCAfit \
		--Predictorcols $prscols \
		--IncludeNofPCs 5)
		
	
	printf "##########\nRun the PRS script (relative to control)\n##########\n"	
		
	PRSsimple_control=$(Rscript	${scriptpath}/prs_compare.R \
		--resultsfile $LCAfit \
		--Predictorcols $prscols \
		--IncludeNofPCs 0 \
		--controlData $controlData)
		
	printf "##########\nRun the PRS script with covariates (relative to control)\n##########\n"	
		
	PRScovar_control=$(Rscript	${scriptpath}/prs_compare.R \
		--resultsfile $LCAfit \
		--Predictorcols $prscols \
		--IncludeNofPCs 5 \
		--controlData $controlData)
		
	printf "##########\nVisualise PRS result\n##########\n"		
		
	#Models must be passed in this order; see script for details
	Rscript	${scriptpath}/extract_prs_model.R \
	$PRSsimple \
	$PRScovar \
	$PRSsimple_control \
	$PRScovar_control
		
end=`date +%s`

runtime=$((end-start))

echo ${runtime}
