### Define a helper function function for adjusting the starts and stscale settings of of a Mplus LCA model


rerunMplusAdjFit<- function(obj,adjustStarts=NULL,adjustStscale=NULL,out=out,
                            datafile=datafile, strts=strts,stscale=stscale){
  
  #Update the starts value in the object with the multiplier defined in adjust starts
  if(!is.null(adjustStarts)){
    obj$ANALYSIS <- gsub(paste0(strts, collapse=" "),paste0(strts*adjustStarts, collapse=" "), obj$ANALYSIS)
    strts <- strts*adjustStarts
  }
  
  #Update the stscale value as + the number specified in adjustStscale
  if(!is.null(adjustStscale)){
    obj$ANALYSIS <- gsub(paste0("STSCALE = ",stscale),paste0("STSCALE = ",stscale+adjustStscale),obj$ANALYSIS)
    stscale <- stscale+adjustStscale
  }

  #Fit model 
  remod <- mplusModeler(obj, modelout=out,
                            dataout = datafile,
                            writeData = 'never',
                            hashfilename = TRUE,
                            run=1L)
  try(print(paste0(out,": ", remod$results$summaries$LL,", starts = ", remod$results$input$analysis$starts,", stscale = ", remod$results$input$analysis$stscale )))
  
  
  #If a saddle point has been reached and the estimator is still the default (MLR) try running the model with MLF estimator
  if(length(grep("THE MODEL ESTIMATION HAS REACHED A SADDLE POINT",remod$results$output))>0 
     && length(grep("ESTIMATOR = MLF",obj$ANALYSIS))==0 ) {
    obj$ANALYSIS <- paste0(obj$ANALYSIS,"\nESTIMATOR = MLF;")
  }
  
  new_replicated <- length(grep("BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED",remod$results$output))>0
  
  
  return(list(
    remod=remod,
    obj = obj,
    new_replicated = new_replicated,
    strts=strts,
    stscale=stscale
  ))
}
