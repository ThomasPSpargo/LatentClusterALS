#### Run across mplus input files to identify the optimum fit for each configuration of data and number of classes

######
## SETUP
######
library(MplusAutomation)

#Read in custom functions from correct directory
helperFuncs <- "/scratch/users/k1802739/ALS_LCA/scripts/mplusHelper"
invisible(lapply(list.files(helperFuncs,full.names = TRUE,pattern=".R"),source))

#Define path to the latent class fits to iterate over (If multiple paths, run for all of them)
runpath<- commandArgs(trailingOnly = TRUE)
 
#Function to loop across paths
mplusPaths <- function(runpath){
  
  #First run initial models, only do so for models without corresponding outfiles
  runModels(runpath,
            recursive=FALSE,replaceOutfile = 'never',
            logFile=paste0(runpath,"/inpmodels_",Sys.Date(),".log"))
  
  models <- readModels(runpath,recursive=FALSE)     #Read in the results of initial model runs
  fewermodels<- dropFinishedModels(models,runpath)  #Use custom function to drop any models which already have final outputs in 'bestout' subdirectory
  
  if(length(fewermodels)>0){
    
    #Identify which vectors have more than 5 classes and may need more starts to find a fit
    lclasses <- vector(mode="numeric")
    for(i in 1:length(fewermodels)){
      size <- fewermodels[[i]]$summaries$NLatentClasses
      if(is.null(size) || size >= 5){lclasses[length(lclasses)+1] <- i}
    }
    
    sclasses <- which(!(fewermodels %in% fewermodels[lclasses])) #identify the smaller classes
    
    #Run this for small model sizes
    lapply(fewermodels[sclasses],findBestMixture,runpath=runpath)
    #Run this with custom starts and high maxstarts limit
    lapply(fewermodels[lclasses],findBestMixture,runpath=runpath,strts = c(1600,400), maxstarts = 15000)
  }
}

#Iterate over runpath values, which indicate the directories to run
lapply(runpath,mplusPaths)
