start.time <- Sys.time()

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(caret)
  library(xgboost)
  library(optparse)
  library(doMC)
  library(MLmetrics) #enables multiclasssummary
})

option_list = list(
  make_option("--rdsDataset", action="store", default=NA, type='character',
              help="Specify the source dataset for which a model will be derived"),
  make_option("--tuneFor", action="store", default=NA, type='character',
              help="Specify the parameter to tune for. Passed to job script"),
  make_option("--assignSuffix", action="store", default="", type='character',
              help="Defines a suffix for the current job; useful if trying different tuning parameters"),
   make_option("--scriptDir", action="store", default=NA, type='character',
               help="Define path to the R script detailing how to process XGBoost fits"),
  make_option("--ncores", action="store", default=30, type='numeric',
              help="Indicate the number of cores that are available for use"),
  make_option("--tuningLogFile", action="store", default="tuning.log", type='character',
              help="Indicate path to a logfile which will be used to store tuning information")
)

opt = parse_args(OptionParser(option_list=option_list))

cat('Options are:\n')
print(opt)

######
### Setup
######

#Obtain directory tree minus model name. This is where output files will be saved
path <- file.path(gsub("_dataset.Rds$","",opt$rdsDataset),paste0("XGBoost_",opt$tuneFor,opt$assignSuffix))
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Directory for the results files
resultspath <- file.path(path,"results")
if(!dir.exists(resultspath)){dir.create(resultspath,recursive = TRUE)}

#Optionally disable writing of tuning via sink to an external log file. 
if(!is.na(opt$tuningLogFile)){
  tuningLog<- file.path(resultspath,opt$tuningLogFile)
  file.create(tuningLog) #Create the log file
  
  cat("Writing tuning details to:", tuningLog,"\n")
} else{
  tuningLog <- NA
}

## read in the nrounds tuning function: tuneNRounds
source(file.path(opt$scriptDir,"xgb_tuneNRounds_func.R"))
## read in the wrapper functions for automated tuning of xgboost in caret
source(file.path(opt$scriptDir,"recursiveTuning_func.R"))

## Read-in data
data<- readRDS(opt$rdsDataset)

## Since classes are unbalanced, determine weighting relative to the class size
classWeights<- min(table(data$C))/table(data$C)

#Loop across records by class weight and assign weighting respective to the class of each individual
rowwiseWeights<- numeric(nrow(data))
for(i in 1:length(classWeights)){
  rowwiseWeights[data$C==names(classWeights)[i]] <- classWeights[i]
}

######
### Configure XGBoost settings
###### 

## Configure general caret training parameters, leaving out 10% each time
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats= 10,
                     search = 'grid',
                     classProbs = TRUE, #Allows inclusion of AUROCC in the curve
                     savePredictions = "final"
) 

#Define the summary function according to the number of classes to model [note that only multi-class is implemented downstream]
if(length(levels(data$C))==2){
  ctrl$summaryFunction <-  twoClassSummary
  objective <- "binary:logistic"
  if(opt$tuneFor=="AUC"){
    opt$tuneFor <- "ROC"
    cat("Performing two class summary - switching caret 'AUC' objective to 'ROC'.\n")
  }
  
} else {
  ctrl$summaryFunction <-  multiClassSummary
  objective <- "multi:softprob"
}

#TUNING PARAMS
# nrounds, ## Boosting Iterations
# max_depth, #Max Tree Depth
# eta, #Shrinkage / learning rate
# gamma, #Minimum Loss Reduction
# colsample_bytree, #Subsample Ratio of Columns #Must be >0
# min_child_weight, #Minimum Sum of Instance Weight)
# subsample #Subsample Percentage #Must be >0

## Create list defining boundaries for possible values of each tuning param, and the steps for values to try in automated retuning of edge params
#Where the limit is infinity, a dummy value of 9999 is set
gridSettings<- list(
  max_depth=c(min=1,max=9999,step=1),
  min_child_weight=c(min=1,max=9999,step=1),
  gamma=c(min=0,max=9999,step=0.1),
  subsample=c(min=0.1,max=1,step=0.05),
  colsample_bytree=c(min=0,max=1,step=0.05), 
  eta=c(min=0,max=1,step=0.02),     #Eta and nrounds do not currently pass through the recursive tuning functions
  nrounds=c(min=1,max=9999,step=10) #Tuning for these is performed in a 2 stage process using xgb.cv functions
)

######
### Perform an initial tuning of nrounds with default settings and a high learning rate (0.3) prior to setting main grid-search params
######

tuneGrid <- data.frame(
  max_depth = 6,
  min_child_weight= 1,
  gamma = 0.0,
  subsample = 1,
  colsample_bytree= 1,
  eta = 0.3, 
  nrounds=100
)

#Extract resampled crossfolds externally using fixed seed for replicability
set.seed(24920)
initResamps <- createMultiFolds(data$C, k = 10,times=10)
#Run xgb.cv and then set the initial number of nrounds
tryCatch({
  if(!is.na(tuningLog)){sink(tuningLog,append = TRUE)}
  initNroundTune<- tuneNRounds(data,rowwiseWeights,tuneGrid,eta=0.3,nrounds=1000,early_stopping_rounds=50,index=initResamps,objective=objective)#k=10,times=10)
  if(!is.na(tuningLog)){sink()}
  
  #Save the tune result
  #saveRDS(initNroundTune, file = file.path(resultspath,"initialNroundTune.Rds"))
  
  #Extract the 3rd quartile for best cutpoint across resamples [3rd quartile selected since this should be close enough to the optimum across resamples]
  setinitNrounds<- ceiling(summary(initNroundTune$stopRows$iter)[["3rd Qu."]])
},error = function(x){
  
  message("Initial nround tuning failed and proceeding with a nominal 100 nrounds Error details parsed as warning below.")
  warning(x)
  
  #If the function fails, nominally select 100 rows as the initial tune
  setinitNrounds <<- 100
})

######
### First stage of proper tuning to configure max_depth, min_child_weight, gamma, holding the others as constant at default settings
######
tuneGrid <- expand.grid(
  max_depth = seq(5,10,by=1),
  min_child_weight= seq(1,6,by=1),
  gamma = seq(0.0,0.4,by=0.1),
  subsample = 1,
  colsample_bytree= 1,
  eta = 0.3, 
  nrounds=setinitNrounds
)

#Set a consistent resampling index across retunes #Using fixed seed for replicability
set.seed(302340)
ctrl$index <- createMultiFolds(data$C, k = 10,times=10)

#run recursive tuning
if(!is.na(tuningLog)){sink(tuningLog,append = TRUE)}
recursiveFit1 <- recursiveWrapTrain(data,opt,tuneGrid,ctrl,rowwiseWeights,ncores=opt$ncores,gridSettings=gridSettings,updateNcores=TRUE)
if(!is.na(tuningLog)){sink()}
xgbfit <- recursiveFit1$finalTrain

######
### Perform 2nd tuning step, tuning subsample and colsample_bytree
######

#Extract the best tune and set the parameter grids across which to train for the second set of tuning features
tuneGrid<- as.list(xgbfit$bestTune)
tuneGrid$colsample_bytree <- tuneGrid$subsample <- seq(0.6,1,by=0.05)
tuneGrid <- expand.grid(tuneGrid)


#Update to a new resampling index across retunes
set.seed(436089)
ctrl$index <- createMultiFolds(data$C, k = 10,times=10)

#run recursive tuning
if(!is.na(tuningLog)){sink(tuningLog,append = TRUE)}
recursiveFit2 <- recursiveWrapTrain(data,opt,tuneGrid,ctrl,rowwiseWeights,ncores=opt$ncores,gridSettings=gridSettings,updateNcores=TRUE)
if(!is.na(tuningLog)){sink()}
xgbfit <- recursiveFit2$finalTrain

######
### Perform 3rd tuning step, tuning eta and nrounds across 2 stages
######
# Note that sink function to log file is (optionally) called externally here

#Extract the best previous tune
tuneGrid<- xgbfit$bestTune

#Update to new resampling seed, which will be used in both stages of the paired nrounds/eta tuning
set.seed(897897)
ctrl$index <- createMultiFolds(data$C, k = 10,times=10)

### Stage 1, using xgb.cv early-stopping functionality determine starting points for a series of values for eta and nrounds
tryEtas<- c(0.01,seq(0.02,0.3,by=0.02))
names(tryEtas) <-  paste0("eta",tryEtas)
if(!is.na(tuningLog)){sink(tuningLog,append = TRUE)}
  tuneEta<- lapply(tryEtas,tuneNRounds,data=data,rowwiseWeights=rowwiseWeights,tuneGrid=tuneGrid,nrounds=1000,early_stopping_rounds=50,index=ctrl$index,objective=objective)#k=10,times=10)
if(!is.na(tuningLog)){sink()}

### Stage 2, using caret, and the xgb.cv-primed partial grid search, tune and optimise
tuneGrid<- xgbfit$bestTune[rep(1,each=length(tryEtas)),]
tuneGrid$eta <- tryEtas
tuneGrid$nrounds <- ceiling(sapply(tuneEta, function(x)summary(x$stopRows$iter)[["Median"]]))


#Run non-recursive tuning
if(!is.na(tuningLog)){sink(tuningLog,append = TRUE)}
  Fit3<- wrapTrain(data, opt, tuneGrid, ctrl, rowwiseWeights, ncores = opt$ncores)
if(!is.na(tuningLog)){sink()}
xgbfit <- Fit3


#Extract the final model and save to file for post-processing
outputFile <- file.path(resultspath,"XGBfinalTune.Rds")
saveRDS(xgbfit,outputFile)


## Generate a comprehensive grid search summary and save to file
## Notes:
# cross-validation folds were resampled for each of search 1, 2, and 3
# search 3 is a partial search paired eta and nrounds, therefore has had 0 retunes in the recursive tuning function
# Pseudorandomisation will cause small discrepancies across recursive retunes
allXGBfits<- rbind(cbind(SearchStage=1,recursiveFit1$allResults),
                   cbind(SearchStage=2,recursiveFit2$allResults),
                   cbind(SearchStage=3,nLoopedRetunes=0,Fit3$results)
                   )

tuneOutpath<- gsub(".Rds$","",outputFile)
if(!dir.exists(tuneOutpath)){dir.create(tuneOutpath)}

write.table(x=allXGBfits,file=file.path(tuneOutpath,"AllTunes.tsv"),sep="\t",row.names = FALSE)

## Save best param config determined across staged retunes to distinct file, including the best performance measure
bestTune<- left_join(xgbfit$bestTune,xgbfit$results)[c(names(xgbfit$bestTune),opt$tuneFor)]
write.table(x=bestTune,file=file.path(tuneOutpath,"BestTune.tsv"),sep="\t",row.names = FALSE)




## Run the summary extraction script on final model
system(paste("Rscript",file.path(opt$scriptDir,"ML_extract.R"),
              "--finalTune",outputFile,
              "--tuneFor",opt$tuneFor,
              "--scriptDir",opt$scriptDir))


# #####
# # Visualise performance across tunes in the staged grid search
# # Not run and not fully implemented; retained for reference.
# #####
# viewTunes <- list()
# 
# plotRoot<- ggplot()+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90))
# 
# #First stage tunes gamma, max_depth, and min_child_weight by recursive grid search
# viewTunes$Tune1<- plotRoot+
#   geom_line(data=recursiveFit1$allResults,aes(x=gamma,y=!!sym(opt$tuneFor),color=as.factor(nLoopedRetunes)))+
#   facet_grid(rows=vars(max_depth),cols=vars(min_child_weight))+
#   labs(caption="Facet rows = max_depth, cols = min_child weight, colour is recursive grid search loops")+
#   scale_colour_discrete(guide="none")
# 
# #Second stage tunes colsample_bytree and subsample by recursive grid search
# viewTunes$Tune2<- plotRoot+
#   geom_line(data=recursiveFit2$allResults,aes(x=colsample_bytree,y=!!sym(opt$tuneFor),color=as.factor(nLoopedRetunes)))+
#   facet_wrap(~subsample)+
#   labs(caption="Facets represent subsamples, and colour is recursive grid search loops")+
#   scale_colour_discrete(guide="none")
# 
# #Third stage involves paired tuning of eta and nrounds
# mutFit3<- Fit3$results %>%
#   mutate(pairedTune=paste0(eta," (",nrounds,")"))
# 
# viewTunes$Tune3 <- plotRoot+
#   geom_line(data=mutFit3,aes(x=pairedTune,y=!!sym(opt$tuneFor),group=1))+
#   labs(x="eta (nrounds)")


##Information
end.time <- Sys.time()
time.taken <- end.time - start.time
cat('Analysis finished at',as.character(end.time),'\n')
cat('Analysis duration was',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n')