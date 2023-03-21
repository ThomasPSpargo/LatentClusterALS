#!/bin/bash
#SBATCH --job-name=PrepLCAsumstats
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-2:00
#SBATCH --ntasks=8
#SBATCH --output=/users/k1802739/SUMSTATS_LCA/%x.log

#Supply trailing argument TRUE to force rerun of completed steps
override=${1:-FALSE}

#Root directory to store all summary statistic files
rootdir="/users/k1802739/SUMSTATS_LCA"
cd $rootdir

#Assumed:
#Presence of the following directory, containing required scripts:
#${rootdir}/scripts
#Note the log file is returned to ${rootdir}

### ALS subset summary statistics are not preprocessed prior to passing to GenoPredPipe. These were shared by van Rheenen et al.; prior preprocessing of the file shared involved renaming columns only and chromosomal position is not included within the file to allow ID positional matching

#Paths to GenoPredPipe 1000 genomes HapMap3 harmonised reference per-chromosome plink file [not in the main directory]
plinkref=/scratch/users/k1802739/GenoPred/GenoPredPipe/resources/data/1kg/1KGPhase3.w_hm3.chr

#Run the preprocessing scripts for each of the publicly available sumstats, if their output file isnt already detected or override is true
mkdir -p ${rootdir}/SZ/
if [[ ! -f ${rootdir}/SZ/SZ_sumstats_PGC.tsv || $override == TRUE ]]; then
	echo "---- Preprocessing SZ ----"
	bash ${rootdir}/scripts/prep_SZ2022_sumstats.sh ${rootdir}/SZ
fi

mkdir -p ${rootdir}/FTD/
if [[ ! -f ${rootdir}/FTD/FTD_GWAS_META.txt || $override == TRUE ]]; then
	echo "---- Preprocessing FTD ----"
	bash ${rootdir}/scripts/prep_FTD_sumstats.sh ${rootdir}/FTD
fi

mkdir -p ${rootdir}/PD/
if [[ ! -f ${rootdir}/PD/PD_sumstats_Nalls.tsv || $override == TRUE ]]; then
	echo "---- Preprocessing PD ----"
	bash ${rootdir}/scripts/prep_PD_sumstats.sh ${rootdir}/PD
fi

mkdir -p ${rootdir}/AD/
if [[ ! -f ${rootdir}/AD/AD_sumstats_Kunkle.txt || $override == TRUE ]]; then
	echo "---- Preprocessing AD ----"
	bash ${rootdir}/scripts/prep_AD_sumstats.sh ${rootdir}/AD
fi

#Space delimited array of GWAS sumstat files to position match with GenoPredPipe reference sumstats
declare -a sumstats=($rootdir/SZ/SZ_sumstats_PGC.tsv $rootdir/AD/AD_sumstats_Kunkle.txt $rootdir/FTD/FTD_GWAS_META.txt $rootdir/PD/PD_sumstats_Nalls.tsv)

#Note that:
#The ID matching protocol recognises the columns and aims to generate a 'complete' SNP column for the main protocol (replacing missingness in the column, or creating one where it is absent (using CHR and BP):
#'SNP','CHR','BP','A1','A2','BETA''FREQ'
 
#loop through the preprocessed summary statistic files
for sum in ${sumstats[@]}; do
		
	#Determine output directory, based on input name, dropping the file extension
	sout=$(echo $sum | sed 's/\..*//')

	echo "------------------------------------------"
	echo "Summary file is: $sum"
	
	#If the outfile doesnt exist, run sumstat ID matcher
	if [[ ! -f $sout.GPPmatch.gz || $override == TRUE ]]; then 
	
		echo "Results will be returned in directory: $(dirname $sout)"
				
		#Check for missing rsIDs and retrieve from reference
		#[NB: THIS STEP ASSUMES POSITIONAL ALIGNMENT between GWAS and REF; some sanity checks within, where possible]
		Rscript ${rootdir}/scripts/sumstat_IDmatcher.R \
			--sumstats $sum \
			--ref_plink_chr $plinkref \
			--writeConflicts TRUE \
			--output ${sout}.GPPmatch \
			--gz T
			
 	else
	  	echo "Outfile already exists for $(basename $sum), no action taken"
	fi
done