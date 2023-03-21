#!/bin/bash
#SBATCH --job-name=submit_mplus_batch
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-0:20
#SBATCH --ntasks=1
#SBATCH --output=/scratch/users/k1802739/ALS_LCA/logs/%x.log

#Sort date into the variable $date, with the format YYYY_MM_DD
printf -v date '%(%Y_%m_%d)T' -1

#Save the path to the log file
logpath=/scratch/users/k1802739/ALS_LCA/logs
mkdir -p ${logpath}/$date #create dated log directory to be used by jobs which will be submitted

#Submit each job file, checking if the job has already completed
for FILE in /scratch/users/k1802739/ALS_LCA/FinalModels/**/*.sh
do 
	echo "The script to run is $FILE";
	
	paths=$(dirname $FILE) #Identify directory path
	
	final_out=$(basename -a ${paths}/bestout/*.out) #Identify corresponding final files

	#For each target directory element
	for inp in ${paths}/*.inp; do
		
		outfile=$(echo $inp | sed s/.inp/.out/g) #Identify corresponding initial outfile
		
		#If the output file doesnt exist, then proceed 
		if [[ ! -f $outfile ]]; then		
			proceed=$"TRUE"
		else #Next, look to the final outputs directory
		
			#Identify which number of clusters is being checked
			final_out_cnumber=$(basename $inp | gawk 'match($0, /.*(c[0-9+]).*/,a) {print a[1]}')
			
			#Set proceed to true, overwrite in loop if outfile match is found
			proceed=$"TRUE"
			
			#Check for the k-clusters match with one of the final outfiles and break loop if found			
			for outs in ${final_out[@]}; do 
				if [[ $outs =~ "$final_out_cnumber" ]]; then
					proceed=$"FALSE"
					break 	
				fi
			done
		fi
		
		#If criteria for $proceed == FALSE have not been met, submit the job
		if [[ $proceed == "TRUE" ]]; then

			echo "Submit the job"
			sbatch -p cpu $FILE

			break #Break the inp loop if job is submitted
		fi

	done # End inp loop
	
	echo "------Filebreak---------"
		
done