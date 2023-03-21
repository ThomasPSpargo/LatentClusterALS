## Define helper function function to drop models which have already been given a 'bestfit' from mplus automation loop

dropFinishedModels <- function(models,runpath){

  #Define function to set a model to null
  elementDrop <- function(mod,runpath){
    
    #Make outpath a subdirectory of inpath
    out_dir <- file.path(runpath,"bestout")
    
    if(dir.exists(out_dir)){
      
      #Look into the outpath for an mplus file with the same root
      outroot<- file.path(out_dir,gsub(".out$","",mod$summaries$Filename))
      
      #See the filenames in the file path
      outfiles<- grep(outroot, list.files(out_dir,full.names=TRUE),value = TRUE)
      
      outfit<- grep(".out$",outfiles,value=TRUE)
      
      #Check to see the existance of outpath file
      if(length(outfit)>0){
        message(paste0("File exists ", mod$summaries$Filename," will not be remodelled"))
        mod <- NULL
      }
    }
    #Return the model, null or unchanged
      return(mod)
  }
  
  #lapply element dropping function
  nullmodels <-lapply(models,elementDrop,runpath=runpath)
  
  #Cut out set as null from the list
  nullmodels<- nullmodels[!sapply(nullmodels,is.null)]
}
