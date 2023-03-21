start.time <- Sys.time()

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  #library(doMC) #Called downstream if running in parallel
  library(randomForest) 
  library(xgboost)
  library(optparse)
  library(kernelshap)
})

option_list = list(
  make_option("--finalTune", action="store", default=NA, type='character',
              help="Specify the source file from which a model will be obtained"),
  make_option("--datasetDir", action="store", default=NA, type='character',
              help="Specify the directory containing all datasets"),
  make_option("--ntasks", action="store", default=1, type='numeric',
              help="Specify number of nodes available for parallelisation"),
  make_option("--nPredictSample", action="store", default=NA, type='numeric',
              help="Indicate approximately how many records to use for feature importance evaluation"),
  make_option("--supplyWeights", action="store", default=TRUE, type='logical',
              help="Indicate TRUE to weight records of background data. Defaults to TRUE"),
  make_option("--outfile", action="store", default=1, type='character',
              help="Path and prefix for the output file")
)

opt = parse_args(OptionParser(option_list=option_list))

options(echo = FALSE)
cat('Options are:\n')
print(opt)


#### Calculate generalisable variable importance metrics with SHAP scores
############################################################

#Read in all the source datasets
dsets<-list.files(opt$datasetDir,pattern="_dataset.Rds",full.names = TRUE)
names(dsets) <-  basename(dsets)
data<- lapply(dsets,readRDS)

#Read in the model for generating shap values
model<- readRDS(opt$finalTune)

#Declare to use parallelisation if there is more than one task
useParallel <- opt$ntasks>1
if(useParallel){
  library(doMC)
  registerDoMC(cores=opt$ntasks)
}

## grep invert T/F to identify which dataset to use for shap prediction [i.e. the training sample for the ML model]:

# whether the multi-class model was used
check<- length(model$levels)
cat("Numlevels:",check,"\n")
notpattern <- check %in% c(4,5)
poss<- grep("class12",names(data),value=TRUE,invert=notpattern)

# size of the training data
check<- nrow(model$trainingData)
cat("NumRows:",check,"\n")
notpattern <- check<5000
poss <- grep("_pheno_",poss,value=TRUE,invert=notpattern)

#Two patterns remain if its genetic or filtered data
if(length(poss)>1){
  #number of features in model; fewer features indicates the sample-matched clinical dataset
  check<- ncol(model$trainingData)
  cat("NumCols:",check,"\n")
  notpattern <- check>10
  poss<- grep("filter",poss,value=TRUE,invert=notpattern)
}
cat("Matched to:",poss,"\n\n")

#Extract the data matching filter criteria
matchedData <- data[[which(names(data)==poss)]]

#Declare the name for the rds file storing the final output
jobFile<- paste0(opt$outfile,"_",class(model$finalModel),"_",gsub("_dataset.Rds","",poss),".Rds")

#Check for the presence of the summary file. If it exists, end the script
if(file.exists(jobFile)){
  cat("Job already complete. Exiting.")
} else {
  
  # Wrap the prediction function
  if(class(model$finalModel)=="xgb.Booster"){
    pfun <- function(object, newdata) {
      predict(object,  newdata = xgb.DMatrix(newdata),reshape=TRUE)
    }
  } else {
    #For a randomForest object
    pfun <- function(object, newdata) {
      predict(object,  newdata = newdata,type="prob")
    }
  }
  
  ## Using stratified sampling via caret functions:
  
  #Determine the proportion of data required to generate ~500 samples for use as background data
  propBG=500/nrow(matchedData)
  partitionBG<- caret::createDataPartition(matchedData$C,p=propBG)[[1]]
  
  #(Optionally) determine the proportion of data required to generate ~opt$nPredictSample samples for SHAP prediction
  if(!is.na(opt$nPredictSample)){
    propPredict=opt$nPredictSample/nrow(matchedData)
    partitionPred<- caret::createDataPartition(matchedData$C,p=propPredict)[[1]]
  } else {
    #if NA, retain all records
    partitionPred <- 1:nrow(matchedData)
  }
  
  #(Optionally) weight the background data by class weights
  if(opt$supplyWeights){
    ## Determine weighting relative to the class size
    classWeights<- min(table(matchedData$C))/table(matchedData$C)
    
    #Loop across records by class weight and assign weighting respective to the class of each individual
    bgWeights<- numeric(length(partitionBG))
    for(i in 1:length(classWeights)){
      bgWeights[matchedData$C[partitionBG]==names(classWeights)[i]] <- classWeights[i]
    }
    cat("Weighting background data by class proportions.\n")
  } else {
    bgWeights=NULL
  }
  
  system.time(
    kshap<- kernelshap(model$finalModel, 
                       X = as.matrix(matchedData[partitionPred,-1]),
                       bg_X = as.matrix(matchedData[partitionBG,-1]),
                       bg_w = bgWeights,
                       parallel = useParallel,
                       pred_fun=pfun)
  )
  
  #Save the SHAP values into an Rds file
  saveRDS(kshap,jobFile)
  
}

# ##Information
end.time <- Sys.time()
time.taken <- end.time - start.time
cat('\nAnalysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')
