library(tidyverse)
library(caret)
library(randomForest)
library(MLmetrics) #enables multiclass summary
library(doMC) #For parallel processing
library(optparse)

#mhttps://towardsdatascience.com/effective-feature-selection-recursive-feature-elimination-using-r-148ff998e4f7

option_list = list(
  make_option("--rdsDataset", action="store", default=NA, type='character',
              help="Specify the source dataset for which a model will be derived"),
  make_option("--tuneFor", action="store", default=NA, type='character',
              help="Specify the parameter to tune for. Passed to job script"),
  make_option("--assignPrefix", action="store", default="", type='character',
              help="Defines a prefix for the current job; useful if trying different tuning parameters"),
  make_option("--ncores", action="store", default=10, type='numeric',
              help="Indicate the number of cores that are available for use"),
  make_option("--scriptDir", action="store", default=NA, type='character',
              help="Define path to directry containing model extraction scripts and required helper funs")
)

opt = parse_args(OptionParser(option_list=option_list))

# opt <- list(rdsDataset="/Users/tom/Downloads/Dec1stmodels/MLprediction/data_pheno_dataset.Rds",tuneFor="AUC",assignPrefix="test")
cat('Options are:\n')
print(opt)

######
### Prepare file structure
######

#Obtain directory tree minus model name. This is where output files will be saved
path <- file.path(gsub("_dataset.Rds$","",opt$rdsDataset),paste0("RF_",opt$tuneFor,opt$assignPrefix))
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

######
### Read-in and configure data
###### 
data<- readRDS(opt$rdsDataset)

#Weightings for sampling to grow trees (as used in Xgboost)

## Since classes are unbalanced, determine weighting relative to the class size
classWeights<- min(table(data$C))/table(data$C)
#Loop across records by class weight and assign weighting respective to the class of each individual
rowwiseWeights<- numeric(nrow(data))
for(i in 1:length(classWeights)){
  rowwiseWeights[data$C==names(classWeights)[i]] <- classWeights[i]
}

#RF style weightings for classwt argument
classwt <- table(data$C)/length(data$C)

######
### Configure random forest settings
###### 

#Control parameters
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats= 10,
                     search = 'grid',
                     classProbs = TRUE,
                     savePredictions = "final"
) 

#For tuning comparability across loops and , resample the same indices a fixed seed
set.seed(428924)
ctrl$index <- createMultiFolds(data$C, k = 10,times=10)

#Define the summary function according to the number of classes to model
if(length(levels(data$C))==2){
  ctrl$summaryFunction <- twoClassSummary
  if(opt$tuneFor=="AUC"){
    opt$tuneFor <- "ROC"
    cat("Performing two class summary - switching caret 'AUC' objective to 'ROC'.\n")
  }
} else {
  ctrl$summaryFunction <- multiClassSummary
}

######
### Configure grid search (internal to and external from caret)
######

## Define hyper-parameter tuning not included within the caret internal grid search
#Test two different tree sizes, but this param shouldnt need tuning since performance reaches an elbow and and won't meaningfully increase past a point
#https://stats.stackexchange.com/questions/348245/do-we-have-to-tune-the-number-of-trees-in-a-random-forest
ntrees <- c(501,1001) 
nodesize <- c(1:40) 
params <- expand.grid(ntrees = ntrees,
                      nodesize = nodesize)

## Define hyper-parameter tuning that IS included with caret
#mtry is the number of variables sampled as candidates in split.
#Upper bound is number of predictor vars; therefore try 1 to the max number
tuneGrid <- expand.grid(.mtry = c(1 : ncol(data)-1))

## Set up the multicore processing
# Note: if the underlying model also uses foreach, the
# number of cores specified above will double (along with the memory requirements)
registerDoMC(cores = opt$ncores)

#### Example roughly followed
# https://stackoverflow.com/questions/57939453/building-a-randomforest-with-caret

#Loop across the external tuning parameters
external_tune <- vector("list", nrow(params))

for(j in 1:nrow(params)){
  #External parameter tuning:
  nodesize <- params[j,2]
  ntree <- params[j,1]
  
  external_tune[[j]] <- train(C~.,
                              data = data,
                              method = "rf",
                              metric = opt$tuneFor,
                              tuneGrid = tuneGrid, #Internal parameter tuning
                              importance = TRUE,
                              classwt = classwt, #Priors of the classes, weights in the model
                              trControl = ctrl,
                              ntree = ntree, #https://stats.stackexchange.com/questions/348245/do-we-have-to-tune-the-number-of-trees-in-a-random-forest
                              weights = rowwiseWeights, #In random forest, these weights are only used for sampling data to grow each tree (not used in any other calculation)
                              nodesize = nodesize)
}

#Assign useful names to list elements
names(external_tune) <- paste("ntrees:", params$ntrees,
                              "nodesize:", params$nodesize)

#Save all results to Rds file
rdsfile <- file.path(path,"RFresults.Rds")
saveRDS(external_tune, file = rdsfile)

#Save best result to distinct Rds file
bestTune<- lapply(external_tune, function(x) x$results[x$results[opt$tuneFor] == max(x$results[opt$tuneFor]),])
finalTune<- which(sapply(bestTune,function(x)x[,opt$tuneFor])==max(sapply(bestTune,function(x)x[,opt$tuneFor])))
model <- external_tune[[finalTune]]

rdsfile <- file.path(path,"RFfinalTune.Rds")
saveRDS(model, file = rdsfile)


## Summarise all tunes in a single table
summaryCols<- c("ntree","nodesize",colnames(external_tune[[1]]$results))
allRFfits = data.frame(matrix(nrow = 0, ncol = length(summaryCols)))
colnames(allRFfits) <- summaryCols

for(i in 1:length(external_tune)){
  allRFfits <- rbind(allRFfits,
                     cbind(ntree=external_tune[[i]]$dots$ntree,
                           nodesize=external_tune[[i]]$dots$nodesize,
                           external_tune[[i]]$results))
}

tuneOutpath<- gsub(".Rds$","",rdsfile)
if(!dir.exists(tuneOutpath)){dir.create(tuneOutpath)}

#Write all tunes to file
allRFfits %>%
  dplyr::arrange(desc(.[,opt$tuneFor])) %>%
  write.table(x=.,file=file.path(tuneOutpath,"AllTunes.tsv"),sep="\t",row.names = FALSE)

#Save best param config to distinct file
bestTune<- allRFfits[allRFfits[,opt$tuneFor]==max(allRFfits[,opt$tuneFor]),c("ntree","nodesize","mtry",opt$tuneFor)]
write.table(x=bestTune,file=file.path(tuneOutpath,"BestTune.tsv"),sep="\t",row.names = FALSE)


##Run the script to extract the best fit
system(paste("Rscript", file.path(opt$scriptDir,"ML_extract.R"),
             "--finalTune",rdsfile,
             "--tuneFor",opt$tuneFor,
             "--scriptDir",opt$scriptDir))

## Visualise performance across tunes
plotTunes<- ggplot(allRFfits,aes(x=nodesize,y=!!sym(opt$tuneFor),color=as.factor(ntree)))+
  geom_line()+
  facet_wrap(~mtry)+
  labs(color="ntree")+
  theme_bw()

ggsave(plot=plotTunes, filename = file.path(tuneOutpath,"RF_tuneChange.pdf"),units="mm",width=150,height=150)

