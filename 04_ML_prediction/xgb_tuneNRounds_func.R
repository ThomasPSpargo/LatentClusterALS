#####
## Helper function for optimising the number of rounds in xgboost with early stopping rounds functionality
## while operating within a caret-style repeatedcv analysis framework
## Please note that I have only accounted for multiclass classification problems (multi:softprob objective), with xgboost maximising AUC
## The the script also assumes that the outcome variable is in the column "C", and that all other columns are features for use in prediction
#####

#### Example application of this function:
# Apply for a single value of eta
#initNroundTune<- tuneNRounds(data,rowwiseWeights,tuneGrid,eta=0.3,nrounds=1000,early_stopping_rounds=50,k=10,times=10)

# For tuning the best eta within lapply: 
# etaTune <- c(0.2,0.3)
#NroundTune_acrossEta<- lapply(etaTune,tuneNRounds,data=data,rowwiseWeights=rowwiseWeights,tuneGrid=tuneGrid,nrounds=1000,early_stopping_rounds=50,k=10,times=10)

## Output:
# stop_summary: A quick measure of the distribution of stopping values recommended by xgboost,
# FullTuneResults: a list of all xgboost result runs
# nroundPlot: (optionally) a ggplot visualising the performance improvements across all resampled cross-fold validations
# stopRows: A data.table summary of all the resamples, iteration stopping points, and auc for the training and test samples at the stopping point
# evalLog: A concatenated summary of train and test xgb.cv evaluation logs across resamples; this used as input data for the ggplot.


# Dependencies
#library(tidyverse)
#library(xgboost)

#### Caret itself doesn't conveniently handle early_stopping_rounds functionality
## #c.f. https://github.com/topepo/caret/issues/641

## Passing directly to Xgboost and using the early stopping functionality identify the optimum Nrounds configuration at a given learning rate
#data - the data to be cross-validated
#rowwiseWeights - the data to be cross-validated
#tuneGrid - the current tuning grid settings used in caret; expects everything to be held constant and nrounds and eta, the parameters tuned here, to be supplied as individual arguments
#eta - the learning rate at which to optimise nrounds; if not specified, the function expects this in the tuneGrid
#nrounds - max number of rounds to try; this is intentionally high by default to encourage early stopping
#early_stopping_rounds - as in xgboost; declare the number of rounds in watchlist dataset without improvement before halting process. If NULL then early stopping will be disabled
#k - as in caret trainControl; number of cross-validation folds
#times - as in caret trainControl; number of repeated cross validations to apply
#incPlot - defaults to FALSE, but set TRUE to return a ggplot of the performance in cross-validated test and train samples across resamples and number of rounds [included as an option and FALSE by default to make ggplot dependency optional]
tuneNRounds <- function(data,rowwiseWeights,tuneGrid,eta=NA,nrounds=1000,early_stopping_rounds=50,k=10,times=10,incPlot=FALSE){
  
  cat("---------\n")
  if(nrow(tuneGrid)>1){
    tuneGrid<- tuneGrid[1,]
    cat("Tuning parameters from tuneGrid argument should only have one row. Subsetting to first row only.\n")
  }
  
  #Assign eta to match the value passed to this function
  if(!is.na(eta)){
    tuneGrid$eta <- eta
  }
  
  cat("Determining the optimum nrounds for learning rate (eta): ", eta,".\nUsing ",k,"-fold cross-validation with ",times," repeated samples.",sep="")
  
  #Always set Nrounds according to the value passed to the Nrounds argument
  tuneGrid$nrounds <- nrounds
  
  #Create repeated cross-validation folds as in caret; the function aims to preserve class balance in the folds
  folds <- createMultiFolds(data$C, k = k,times=times)
  
  #Invert the folds to identify the test data
  testFolds<-lapply(folds,function(x,y) y[!y %in% x],1:nrow(data))
  #intersect(testFolds[[1]],folds[[1]]) #Intersect to check that there is no overlap
  
  #Generate a parameter list for direct XGBoost tuning
  params<- c(tuneGrid[names(tuneGrid)!="nrounds"],
             objective = "multi:softprob",
             eval_metric = "auc",
             num_class=length(levels(data$C)),
             nthread=1)
  
  #Internal function for formatting datasets, distinguishing the clusters 'C' from the predictors, and assigning weights
  formatData<- function(x,weights){
    y <- xgb.DMatrix(as.matrix(x[names(x)!="C"]),
                     label=as.numeric(x$C)-1, #xgboost: Class is represented by a number and should be from 0 to num_class - 1.
                     weight=weights
                     )
    return(y)
  }
  xgbdata<- formatData(data,rowwiseWeights)
  
  #Loop across resamples and apply xgb.cv across folds
  resamples<-unique(gsub("Fold[0-9]+\\.","",names(testFolds)))
  resampXGB_CV<-vector(mode="list",length=times)
  names(resampXGB_CV) <- resamples
  for(i in 1:length(resamples)){
    tFlds<- testFolds[grep(resamples[i],names(testFolds))]
    
    #### xgb.cv should only ger partitions for each resample.
    resampXGB_CV[[i]]<- xgb.cv(
      data = xgbdata,
      params = params,
      nrounds=tuneGrid$nrounds,
      verbose = 0,
      folds=tFlds,
      prediction=TRUE,
      early_stopping_rounds = early_stopping_rounds,
      maximize = TRUE,
    )
  }
  
  #Extract the stopping thresholds
  stopPoints<- sapply(resampXGB_CV,function(x){x$best_ntreelimit})
  
  #Summarise the distribution of stop points
  stop_summary<- summary(stopPoints)
  
  cat("Summary for best number of rounds indicated by xgb.cv across repeated resamples:\n")
  print(stop_summary)
  cat("---------\n")
  
  #Extract the rows selected as stop points and convert to a tabular summary
  stopRows<- mapply(function(x,row){x$evaluation_log[row,]},
                    resampXGB_CV,stopPoints,SIMPLIFY=FALSE)
  stopRows <- cbind(resample=names(stopRows),do.call(rbind, stopRows))
  
  #Extract the nrounds performance by resample, for plotting.
  nroundPerform <- data.frame(resample=character(0),iter=numeric(0),train_auc=numeric(0),test_auc=numeric(0))
  for(i in 1:length(resampXGB_CV)){
    
    #Extract the resample ID and the evaluation log, combined into a single dataframe
    indivcurve<- cbind(resample=rep(names(resampXGB_CV)[i],times=nrow(resampXGB_CV[[i]]$evaluation_log)),
                       resampXGB_CV[[i]]$evaluation_log)
    
    nroundPerform<- rbind(nroundPerform,
                          indivcurve)
    
  }
  
  if(incPlot){
    nroundPlot<- nroundPerform %>%
      pivot_longer(.,cols=c("test_auc_mean","train_auc_mean"),values_to = "AUC",names_to = "dataset") %>%
      ggplot(.,aes(x=iter,y=1-AUC,alpha=resample,color=dataset))+
      geom_line(linewidth=1.5)+
      scale_color_manual(values=c("test_auc_mean"="red","train_auc_mean"="blue"),labels=c("Test","Train"))+
      theme_bw()+
      geom_vline(data=as.data.frame(t(stop_summary)),aes(xintercept=Freq,lty=Var2),alpha=0.6)+
      labs(lty = "Overall\nstop cutoffs",color="Resampled\ndataset\npartition",alpha="",x="Number of rounds",title=paste0("XGBoost across resamples with eta = ",unique(tuneGrid$eta)))
    
    
  } else {
    nroundPlot <- "Set the incPlot=TRUE to return a ggplot summarising performance across nrounds"
  }
  
  #Return details of the tuning results
  return(list(
    stop_summary=stop_summary,
    FullTuneResults=resampXGB_CV,
    nroundPlot=nroundPlot,
    stopRows=stopRows,
    evalLog=nroundPerform
  ))
  
} #End Nround tuning function