### Specify a function for recursively running mplus with tweaked settings in order to find the best-fit latent class model for a given dataset and number of classes.


suppressMessages(library(MplusAutomation))

#Users can optionally specify the starting values for starts, stscale and maxstarts
#Maxstarts is the value used as a breakpoint when iterating to try and find a good fit
findBestMixture <- function(mod,runpath,
                            strts=c(400,100),
                            stscale=5,
                            maxstarts=5000) {
  
  ##### Initial prep
  old <- Sys.time()       #get start time
  startdate <- Sys.Date() #get start date
  
  #Find the full file-path for the input file and declare as the current best model
  out_best <- inpath <- file.path(runpath,mod$summaries$Filename)
  
  out_dir <- file.path(runpath,"bestout") #Make outpath a subdirectory of runpath
  if(!dir.exists(out_dir)){dir.create(out_dir,recursive = TRUE)}
  
  mid_dir <- file.path(runpath,"allfits") #Make middle-step directory a subdirectory of runpath
  if(!dir.exists(mid_dir)){dir.create(mid_dir,recursive = TRUE)}
  
  #Name columns of summary file
  sumcols <- matrix(c("title",
                      "estimator",
                      "startvalues",
                      "stscale",
                      "bestll_allruns",
                      "prevbestLL",
                      "finalLL",
                      "AIC",
                      "BIC",
                      "aBIC",
                      "entropy",
                      "class_counts",
                      "smallestCProb",
                      "prevrep",
                      "newrep",
                      "warnings",
                      "errors",
                      "elapsed",
                      "sanity_check"
  ),nrow=1)
  
  #If no summary file currently exists, create one with the required headers
  if(!file.exists(file.path(out_dir,paste0("bestfits",startdate,".csv")))){
    
    write.table(sumcols, 
                file=file.path(out_dir,paste0("bestfits",startdate,".csv")),
                sep=",",
                row.names = FALSE,
                col.names = FALSE)
    
  }
  
  if(is.null(mod$summaries$NLatentClasses)){
    #Sanity check step
    warning("Model has no latent classes")
    
  } else if(mod$summaries$NLatentClasses==1){
    #In 1 class models the solution is found immediately, copy the original file into end directory
    message("Model has 1 class")
    
    #Copy the files into the 'bestout' summary directory
    file.copy(c(inpath, gsub(".out$",".inp",inpath)),out_dir)
    
    #Create across 1-class model summary 
    #Define empty matrix
    prgrss <- matrix(rep(NA_character_,length(sumcols)),ncol=length(sumcols))
    
    #Insert values into the correct element of the prgrss matrix
    try(prgrss[1] <-  mod$input$title)
    #try(prgrss[2] <-  mod2$results$summaries$Estimator)
    #try(prgrss[3] <-  paste(strts, collapse="->"))
    #try(prgrss[4] <-  stscale)
    #try(prgrss[5] <- bestll)
    #try(prgrss[6] <-  prev)
    try(prgrss[7] <-  mod$summaries$LL) 
    try(prgrss[8] <-  mod$summaries$AIC)
    try(prgrss[9] <-  mod$summaries$BIC)
    try(prgrss[10] <-  mod$summaries$aBIC)
    try(prgrss[11] <- 1) #Entropy must be 1 because only one class
    try(prgrss[12] <- mod$class_counts$mostLikely$count)
    #try(prgrss[13] <- min(diag(mod$class_counts$avgProbs.mostLikely))) #Smallest class probability
    #try(prgrss[14] <- prev_replicated)
    #try(prgrss[15] <- new_replicated)
    try(prgrss[16] <- paste(sapply(mod$warnings, paste, collapse=" "), collapse=" ->"))
    try(prgrss[17] <- paste(sapply(mod$errors, paste, collapse=" "), collapse=" ->"))
    #try(prgrss[18] <- time)
    #try(prgrss[19] <- sanitycheck,silent = TRUE)
    
    write.table(prgrss, 
                file=file.path(out_dir,paste0("bestfits",startdate,".csv")),
                sep=",",
                row.names = FALSE,
                col.names = FALSE,
                append=TRUE
    )
    
  } else { #i.e. if(mod$summaries$NLatentClasses>1)
    
    
    ##########################################
    # extract relevant info from initial fit #
    ##########################################           
    sanitycheck <- NULL #Set the sanity check parameter to NULL
    
    #Define current best ll value (to be added to the output) - conditional upon the presence of a successful fit in the data
    if(!is.null(mod$summaries$LL)){bestll <- mod$summaries$LL   
    } else {bestll <- NULL}
    prev <- bestll #Take the loglikelihood of this model as a reference point for indication of fit
    
    prev_replicated <- length(grep("BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED",mod$output))>0 #Check whether the new best value is replicated
    
    #Create mod best from the first model attempt
    mod_best <- list(results = mod,
                     ANALYSIS = paste0("STARTS =", mod$input$analysis$starts,";\nSTSCALE = ",mod$input$analysis$stscale,";")
    )
    
    
    #Based on input path, create new filepath to be used for the intermediate outputs, where runs will be tested iteratively:
    out <- sub(".out", "_re.inp",file.path(mid_dir,mod$summaries$Filename)) #sub is used to convert to inp file string, for use in MplusModeller
    
    #Create a general purpose datafile directory for use in repeated loops
    dataout<- file.path(mid_dir, gsub(".csv",".dat",basename(mod$input$data$file)))
    
    #####
    ###  Prepare mplus object for reruns
    #####
    
    
    #Convert the Mplus model into a mplusAutomation mplusObject and write the datafile to the output directory
    prep <- prepObject(mod,out,dataout)
    obj <- prep$obj
    datafile <- prep$dataout
    
    #Update the starts and scale parameters according to value used in the rerun
    obj$ANALYSIS <- sub(paste0("STARTS = ", mod$input$analysis$starts, ";"),
                        paste0("STARTS = ", paste0(strts, collapse=" "), ";"), obj$ANALYSIS)
    
    obj$ANALYSIS <- sub(paste0("STSCALE = ",mod$input$analysis$stscale,";"),
                        paste0("STSCALE = ",stscale,";"), obj$ANALYSIS)
    
    #Run the analysis
    rerun<- rerunMplusAdjFit(obj,adjustStarts=NULL,adjustStscale=NULL,out=out,
                             datafile=datafile, strts=strts,stscale=stscale)
    
    #Extract elements to be checked from output list
    mod2<- rerun$remod
    new_replicated <- rerun$new_replicated
    obj <- rerun$obj
    strts <- rerun$strts #unchanged here
    stscale <- rerun$stscale #unchanged here
    
    
    #If a new best ll is found, replace bestll with this value, and keep the best model details
    if((!is.null(mod2$results$summaries$LL) && is.null(prev)) ||
       (!is.null(mod2$results$summaries$LL) && bestll <= mod2$results$summaries$LL)){
      bestll <- mod2$results$summaries$LL
      mod_best <- mod2
      out_best <- out
    } 
    
    
    #Sanity check to test result using a wider stscale (i.e. reduce chance of getting stuck with local maxima)
    #================================
    #Rename the mplus .inp and .out files by number of loops,
    out_wide <- file.path(dirname(out),gsub("re.*",paste0("re_stscale",stscale+10,".inp"), basename(out)))
    
    #Rerun the analysis
    rerun<- rerunMplusAdjFit(obj,adjustStarts=NULL,adjustStscale=10,out=out_wide,
                             datafile=datafile, strts=strts,stscale=stscale)
    
    #Extract elements to be checked from output list
    mod_wide <- rerun$remod
    obj_wide <- rerun$obj
    #new_replicated <- rerun$new_replicated
    #strts <- rerun$strts
    #stscale <- rerun$stscale
    
    #If the fit is better or the same, or if the original model is yielding null
    if((!is.null(mod_wide$results$summaries$LL) && mod2$results$summaries$LL <= mod_wide$results$summaries$LL)
       || is.null(mod2$results$summaries$LL)){
      
      out <- out_wide
      
      stscale <- rerun$stscale #Update stscale to reflect that widening appears improved the fit
      
      #Save new best loglikelihood as new and old wide, then try additional widths
      new_wide <- mod_wide$results$summaries$LL
      old_wide <- mod_wide$results$summaries$LL
      
      mod_wide_best <- mod_wide #Indicate this is the best 'wider' fit
      
      #Object to test whether new-wide repeatedly does not improve (trigger break when dup>=2)
      dup <- 0 
      
      #So long as new_wide is improving the solution, loop models
      while(old_wide <= new_wide){
        
        out_wide <- file.path(dirname(out),gsub("re.*",paste0("re",strts[1],"stscale_",stscale+10,".inp"), basename(out))) #Rename the mplus .inp and .out files by width of scale
        
        #Rerun the analysis
        rerun<- rerunMplusAdjFit(obj=obj_wide,adjustStarts=NULL,adjustStscale=10,out=out_wide,
                                 datafile=datafile, strts=strts,stscale=stscale)
        
        #Extract elements to be checked from output list
        mod_wide <- rerun$remod
        obj_wide <- rerun$obj
        
        #If the same or better, update the stscale and update the loop
        if(!is.null(new_wide) && old_wide <= new_wide) stscale <- rerun$stscale
        
        #If the new value actually exceeds the previous, overwrite the out value
        if(!is.null(new_wide) && old_wide < new_wide){
          old_wide <- new_wide
          mod_wide_best <- mod_wide #Save the overall model
          dup <- 0 #Reset dup if a new best value is hit
        }
        
        if(!is.null(new_wide) && old_wide == new_wide)dup <- dup+1 #If new-wide yields repeated results +1 to dup, break loop once dup hits +2
        
        if (stscale >50 || is.null(new_wide) || dup>=2) {
          print("stscale loop break")
          break
        }
      }
      
      #Once loop breaks, take the best solution
      mod_best <- mod2 <- mod_wide_best
      
      obj_wide$ANALYSIS <- mod2$ANALYSIS #Match STSCALE to the one from the best model
      
      stscale <- as.numeric(mod2$results$input$analysis$stscale) #Reset this correctly
      obj <- obj_wide 
      out_best <- out 
      bestll <- mod2$results$summaries$LL #replace bestll and update best models
      
      #Extract the fit summaries for comparison to the subsequent 2*starts fit
      prev <- mod2$results$summaries$LL #Extract the new loglikelihood comparison value
      prev_replicated <- length(grep("BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED",mod2$results$output))>0 #Check whether the new best value is replicated
      
      #Rename the mplus .inp and .out files
      out <- file.path(dirname(out),gsub("re.*",paste0("re",strts[1]*2,"_stscale_",stscale,"_trydbl.inp"),basename(out)))
      
      #Rerun the analysis with double the number of starts
      rerun<- rerunMplusAdjFit(obj,adjustStarts=2,adjustStscale=NULL,out=out,
                               datafile=datafile, strts=strts,stscale=stscale)
      
      #Extract elements to be checked from output list
      mod2 <- rerun$remod
      new_replicated <- rerun$new_replicated
      obj <- rerun$obj
      strts <- rerun$strts
      #stscale <- rerun$stscale
      
      #If a new best ll is found, replace bestll with this value, and keep the best model details
      if(!is.null(mod2$results$summaries$LL) && bestll < mod2$results$summaries$LL){
        bestll <- mod2$results$summaries$LL
        mod_best <- mod2
        out_best <- out
      }
      
    } #End check of wider stscales
    
    
    #From here, the script iterates on the  mod2 starts values to try and achieve ML    
    #=================================
    #If the matching criteria are met and/or the sanity check not passed,
    #begin a while loop until break point is met (initial starts > 5000)
    
    while(!(identical(mod2$results$summaries$LL, prev) 
            && !is.null(mod2$results$summaries$LL) #Additional check to make sure null isn't reached
            && new_replicated
            && prev_replicated
            #            && length(grep("FIXED TO AVOID SINGULARITY",mod2$results$output))==0
    )){
      
      #Extract elements for comparison from previous run:
      prev_replicated <- new_replicated
      prev <- mod2$results$summaries$LL
      
      #Rename the mplus .inp and .out files.
      out <- file.path(dirname(out),gsub("re.*",paste0("re",strts[1]*2,".inp"), basename(out)))
      
      #Rerun the analysis with double the number of starts
      rerun<- rerunMplusAdjFit(obj,adjustStarts=2,adjustStscale=NULL,out=out,
                               datafile=datafile, strts=strts,stscale=stscale)
      
      #Extract elements to be checked from output list
      mod2 <- rerun$remod
      new_replicated <- rerun$new_replicated
      obj <- rerun$obj
      strts <- rerun$strts
      #stscale <- rerun$stscale
      
      #If a new best ll is found, replace bestll with this value, and keep the best model details
      if(!is.null(mod2$results$summaries$LL) && bestll < mod2$results$summaries$LL){
        bestll <- mod2$results$summaries$LL
        mod_best <- mod2
        out_best <- out
      }
      
      #Stop the loop if starts exceeds the threshold of 5000 starts
      #Also stop the loop if both of the previous two models are returning at NA, therefore try progressing onto the next step
      if (strts[1] >= maxstarts || is.null(mod2$results$summaries$LL) && is.null(prev)) {
        break
      }
    }
    
    
    #From here, the script prepared output files
    #=================================
    #calculate elapsed time
    time <- Sys.time() - old # calculate difference
    
    #If the best loglikelihood is a better fit than the final model- then store the best record (in addition to the final record)
    if(!is.null(mod2$results$summaries$LL) && !is.null(bestll) && bestll > mod2$results$summaries$LL){
      #Find and match the pattern for number of starts
      try(m<- regexpr("STARTS.*0;", mod_best$ANALYSIS))
      try(pat_starts <- regmatches(mod_best$ANALYSIS,m))
      #Find and match pattern for stscale
      try(m<- regexpr("STSCALE.*;", mod_best$ANALYSIS))
      try(pat_scale <- regmatches(mod_best$ANALYSIS,m))
      
      #Create empty matrix for summary of results
      prgrss <- matrix(rep(NA_character_,length(sumcols)),ncol=length(sumcols))
      
      #Insert values into the correct element of the progress matrix
      try(prgrss[1] <-  mod$input$title)
      try(prgrss[2] <-  mod_best$results$summaries$Estimator)
      try(prgrss[3] <-  pat_starts)
      try(prgrss[4] <-  pat_scale)
      try(prgrss[5] <- bestll)
      ##try(prgrss[6] <-  prev)
      try(prgrss[7] <-  mod_best$results$summaries$LL)
      try(prgrss[8] <-  mod_best$results$summaries$AIC)
      try(prgrss[9] <-  mod_best$results$summaries$BIC)
      try(prgrss[10] <-  mod_best$results$summaries$aBIC)
      try(prgrss[11] <- mod_best$results$summaries$Entropy)
      try(prgrss[12] <- paste0(mod_best$results$class_counts$mostLikely$count,collapse =": "))
      try(prgrss[13] <- min(diag(mod_best$results$class_counts$avgProbs.mostLikely))) #Smallest class probability
      #try(prgrss[14] <- prev_replicated)
      try(prgrss[15] <- length(grep("BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED",mod_best$results$output))>0)
      try(prgrss[16] <- paste(sapply(mod_best$results$warnings, paste, collapse=" "), collapse=" ->"))
      try(prgrss[17] <- paste(sapply(mod_best$results$errors, paste, collapse=" "), collapse=" ->"))
      try(prgrss[18] <- time)
      try(prgrss[19] <- sanitycheck,silent = TRUE)
      
      write.table(prgrss, 
                  file=file.path(out_dir,paste0("bestfits",startdate,".csv")),
                  sep=",",
                  row.names = FALSE,
                  col.names = FALSE,
                  append=TRUE
      )
      
      #Copy the files into the 'bestout' summary directory
      file.copy(c(out_best, gsub(".inp$",".out",out_best)),out_dir)
    }
    
    #Create empty matrix for summary of results
    prgrss <- matrix(rep(NA_character_,length(sumcols)),ncol=length(sumcols))
    
    #Insert values into the correct element of the progress matrix
    try(prgrss[1] <-  mod$input$title)
    try(prgrss[2] <-  mod2$results$summaries$Estimator)
    try(prgrss[3] <-  paste(mod2$results$input$analysis$starts, collapse="->"))
    try(prgrss[4] <-  mod2$results$input$analysis$stscale)
    try(prgrss[5] <-  bestll)
    try(prgrss[6] <-  prev)
    try(prgrss[7] <-  mod2$results$summaries$LL)
    try(prgrss[8] <-  mod2$results$summaries$AIC)
    try(prgrss[9] <-  mod2$results$summaries$BIC)
    try(prgrss[10] <-  mod2$results$summaries$aBIC)
    try(prgrss[11] <- mod2$results$summaries$Entropy)
    try(prgrss[12] <- paste0(mod2$results$class_counts$mostLikely$count,collapse =": "))
    try(prgrss[13] <- min(diag(mod2$results$class_counts$avgProbs.mostLikely))) #Smallest class probability
    try(prgrss[14] <- prev_replicated)
    try(prgrss[15] <- new_replicated)
    try(prgrss[16] <- paste(sapply(mod2$results$warnings, paste, collapse=" "), collapse=" ->"))
    try(prgrss[17] <- paste(sapply(mod2$results$errors, paste, collapse=" "), collapse=" ->"))
    try(prgrss[18] <- time)
    try(prgrss[19] <- sanitycheck,silent = TRUE)
    
    
    message(paste0(prgrss[,c(1:7,14,15,18)],collapse = " "))
    
    write.table(prgrss, 
                file=file.path(out_dir,paste0("bestfits",startdate,".csv")),
                sep=",",
                row.names = FALSE,
                col.names = FALSE,
                append=TRUE
    )
    
    #Copy the files into the 'bestout' summary directory
    file.copy(c(out, gsub(".inp$",".out",out)),out_dir)
    
  } #End check of class sizes
} #end mixtures function 
