#   library(tidyverse)
#   library(optparse)
#   library(caret)
#   library(xgboost)
#   library(optparse)
#   library(doMC)
#   library(MLmetrics) #enables multiclasssummary


#### This script defines two wrapper functions for running caret and xgboost with recursive tuning updates:

#wrapTrain is a wrapper for caret's train function [implemented currently for xgbTree only]
#recursiveWrapTrain is a wrapper for wrapTrain which uses a while loop to automate tuned parameter optimisation if occurring at grid search edges

wrapTrain <- function(data,opt,tuneGrid,ctrl,rowwiseWeights,ncores=1){
  
  #Return an info message describing the grid search parameters which are being tuned and held constant
  gridSize<- sapply(tuneGrid,function(x)length(unique(x)))
  tunedInfo<- tuneGrid[gridSize>1]
  tunedInfo<- lapply(tunedInfo,unique)
  constantInfo<- tuneGrid[gridSize==1]
  
  cat("------------\nTuning parameters with grid search size: ",nrow(tuneGrid),". Trying:\n",paste0(names(tunedInfo)," = [",lapply(tunedInfo,paste,collapse=", "),"]",collapse=",\n"),".\n",sep='')
  cat("\nParams held constant are:\n",paste0(names(constantInfo)," = [",unique(constantInfo),"]",collapse=",\n"),".\n",sep='')
  
  if(ncores>1){
    registerDoMC(cores = ncores)
    ctrl$allowParallel <- TRUE
    cat('Using',ncores,'cores.\n')
  } else {
    ctrl$allowParallel <- FALSE
    cat('Using sequential processing.\n')
  }
  
  start.time <- Sys.time()
  cat("Starting at:",as.character(start.time),"\n")
  
  #Train the model
  #Inherit dataframe and tuning parameters from input. Only necessary column is the factor class 'C'
  fit <- train(C~.,
                  data = data, 
                  method = "xgbTree",
                  metric= opt$tuneFor, #Caret's interpretation of optimisation metric
                  eval_metric = "auc", #xgboost's interpretation of metric
                  nthread = 1, #Force xgboost to use 1 thread and solely use parallel processing implemented by caret (cf. https://stackoverflow.com/questions/39528392/parallel-processing-with-xgboost-and-caret)
                  tuneGrid = tuneGrid,
                  #early_stopping_rounds=50, This feature does not work with the default implementation: https://github.com/topepo/caret/issues/641
                  maximize=TRUE, #Maximise MUST be set alongside early stopping rounds, indicating that higher eval score is better
                  trControl = ctrl,
                  weights= rowwiseWeights, #Caret xgbTree documentation indicates that weights (wts) parameter will be passed to XGB from caret. 
                  verbosity=0 #Silence warning about ntree limit [A depreciated xgboost argument that hasn't yet been implemented in caret]
  )

  end.time <- Sys.time()
  cat("Done! Best params:\n")
  print(fit$bestTune)
  
  cat("Finishing at:",as.character(end.time),"\n")
  time.taken <- end.time - start.time
  cat('Time taken for tuning was ',as.character(round(time.taken,2)),attr(time.taken, 'units'),'\n------------\n')
  
  return(fit)
}


#### Recursive wrapping of wrapTrain to permit automatic retuning at edges
recursiveWrapTrain <- function(data,opt,tuneGrid,ctrl,rowwiseWeights,ncores=1,gridSettings,updateNcores=TRUE,maxloops=5){
  
  #Store all caret train objects in a list
  allTrains <- list()
  
  #Perform the first tune
  fit <- wrapTrain(data,opt,tuneGrid,ctrl,rowwiseWeights,ncores=ncores)
  
  #Save to list
  allTrains[[1]] <- fit
  
  #Save the results for concatenation with any repeated loops
  allResults<- cbind(nLoopedRetunes=0,fit$results)
  
  #Identify the parameters that are (still) being tuned
  tunedParams<- names(tuneGrid[sapply(tuneGrid,function(x)length(unique(x)))>1])
  
  #Counting object to identify current number of retunes (if any)
  nRetunes <- 0
  
  #Check whether any of the tune parameters are at the edges of the grid search, and if they are and are not at the defined min or max, repeat tuning
  while(any(
    mapply(function(x,y){x %in% y}, #Logic check for whether tuned params are at the edge of the tuneGrid
           as.numeric(fit$bestTune[tunedParams]), #X is the best-tune values of the tuneGrid
           lapply(tuneGrid[tunedParams],range) #Y is the min and max values of the tuneGrid
    ) &
    mapply(function(x,y){!x %in% y}, #Logic check for whether tuned params are NOT already at their defined limits
           as.numeric(fit$bestTune[tunedParams]), #X is the best-tune values of the tuneGrid
           lapply(gridSettings[tunedParams],range) #Y is the min and max values allowed in the tuneGrid, as defined in gridSettings
    )
  )){
    
    #Count the number of retunes performed
    nRetunes <- nRetunes+1
    cat("############\n Commencing retune number: ",nRetunes,"\n")
    
    #Create a retuning grid based on the best tune
    reTune<- as.list(fit$bestTune)
    
    #Indicate tune scale of further stepping for any expansions of the tuned parameters
    #Param info: https://xgboost.readthedocs.io/en/latest/parameter.html.
    #Some tuning guidance: https://www.analyticsvidhya.com/blog/2016/03/complete-guide-parameter-tuning-xgboost-with-codes-python/
    tunescale<-sapply(gridSettings[tunedParams],function(x)x[["step"]])
    
    ### Loop across tune params and extend any where the best fit occurred at the edges that are not also parameter limits
    for(i in 1:length(tunescale)){
      if(fit$bestTune[[tunedParams[i]]]==min(tuneGrid[tunedParams[i]]) && fit$bestTune[[tunedParams[i]]]!=gridSettings[[tunedParams[i]]][["min"]]
         ){
        #Try some slightly lower values
        reTune[[tunedParams[i]]] <- fit$bestTune[[tunedParams[i]]]+seq(-3*tunescale[i],tunescale[i],by=tunescale[i])
      } else if (fit$bestTune[[tunedParams[i]]]==max(tuneGrid[tunedParams[i]]) && fit$bestTune[[tunedParams[i]]]!=gridSettings[[tunedParams[i]]][["max"]]
                 ){
        #Try some slightly higher values
        reTune[[tunedParams[i]]] <- fit$bestTune[[tunedParams[i]]]+seq(-tunescale[i],3*tunescale[i],by=tunescale[i])
      }
    }
    
    #Prune the reTune grid to retain only values within the grid limits
    reTune[tunedParams] <- mapply(function(x,y){x[x>=y["min"] & x<=y["max"]]},
                                  reTune[tunedParams], #X is the best-tune values of the tuneGrid
                                  gridSettings[tunedParams], #Y is the min and max values allowed in the tuneGrid, as defined in gridSettings
    SIMPLIFY=FALSE)
    
    #Overwrite the tuneGrid
    tuneGrid <- expand.grid(reTune)
    
    #Update the parameter indicating what is being tuned
    tunedParams<- names(tuneGrid[sapply(tuneGrid,function(x)length(unique(x)))>1])
    
    #End the loop if only an a single gridsearch row is returned
    if(nrow(tuneGrid)==1){
      cat("Terminate while loop early, grid has 1 row.\n")
      break
    }
    
    if(updateNcores){
      #Reassign multicores, using a single core per row of gridsearch if multiple cores are available.
      #Otherwise, calculate the optimum number of cores
      Reset_cores <-nrow(tuneGrid)
      if(Reset_cores>ncores){
        npercore<- Reset_cores/ncores
        Reset_cores<- ceiling(ncores/npercore)
      }
      cat("Registering cores: ",Reset_cores,"\n")
    } else { 
      Reset_cores <- ncores
      cat("Holding ncores constant at: ",Reset_cores,"\n")
    }
    
    #Perform the retune as necessary
    fit <- wrapTrain(data,opt,tuneGrid,ctrl,rowwiseWeights,ncores=Reset_cores)
    
    #Save the train object
    allTrains[[length(allTrains)+1]] <- fit
    
    #Combine the grid search performance with the earlier results
    allResults<- rbind(allResults,
                       cbind(nLoopedRetunes=nRetunes,fit$results)
                       )
    
    if(maxloops==nRetunes){ #Force stop if there have been maxloops retunes - arbitrarily set to the default of 5
      warning(nRetunes,"repeat tunes occurred, stopping the loop.\n")
      break
    } 
  }
  return(list(finalTrain=fit,allTrains=allTrains,allResults=allResults,tuneGrid=tuneGrid))
}