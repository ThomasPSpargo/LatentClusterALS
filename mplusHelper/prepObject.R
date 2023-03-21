### Define a function to prepare an Mplus object to be used in MplusAutomation and perform a dry run of analysis in the output directory

suppressPackageStartupMessages({
  library(MplusAutomation)
  library(data.table)
  library(dplyr)
})

#mod is the MplusAutomation Mplus model object
#out is the file to write .out to
#dataout is the name of the .dat file where the dataset will be saved
prepObject <- function(mod,out,dataout,newData=NULL,hashfilename=TRUE){
  
  if(!is.null(newData)){ dat <- newData #Optionally supply new dataset
  } else {
    
    #obtain the variable names for the model, using custom function
    spreadnames<- spreadMplusHyphens(mod$input$variable$names)
    
    #Identify input file, gsub any blank spaces in directory steps, caused by linebreaks. Use two gsubs to capture spaces before and after '/' separators.
    infile <- gsub("/ +","/",mod$input$data$file) %>%
      gsub(" +/","/",.)
    
    #Import dataout and assign column names
    dat <- data.table::fread(file=infile,
                             header=FALSE,
                             col.names = spreadnames,
                             na.strings = ".")
  } 
  
  usevars <- spreadMplusHyphens(mod$input$variable$usevariables)
  
  mis_col<- usevars[which(!(usevars %in% colnames(dat)))] #Check whether all usevars are present within the cols of dat
  
  if(length(mis_col)>0){ #If any of the usevariables are missing from the target dataset, declare in message and drop them from usevariables
    
    warning("Mplus expects the following input usevariables:", paste(usevars,collapse = " "),".\n",
            "The dataset provided is missing the following:",paste(mis_col,collapse= " "),".\n",
            "usevariables has been adjusted to exclude the missing columns; please check that this is acceptable in your use-case."
    )
    
    usevars <- usevars[which(usevars %in% colnames(dat))] #Keep only values present in the new dataset
  }
  
  
  #In case many variables are included in model, break the string if it hits greater than 70 character length
  usevars<- stringBreak(paste0(usevars,collapse=" "),sep = " ",buffer = 70)
  
  
  if(!is.null(mod$input$variable$useobservations)){useobs<- paste0("USEOBSERVATIONS = ",mod$input$variable$useobservations,";\n")
  } else {useobs <- ""}
  
  #Build mplus object based on mod (a mplus.list class object)
  obj <- mplusObject(
    TITLE= mod$input$title,
    rdata = data.frame(dat),
    VARIABLE =
      paste0("USEVARIABLES ARE ", usevars,";\n", 
             useobs,
             "IDVARIABLE is ",mod$input$variable$idvariable, ";\n",
             "NOMINAL = ",mod$input$variable$nominal, ";\n",
             "CATEGORICAL = ",mod$input$variable$categorical, ";\n",
             "SURVIVAL = ", mod$input$variable$survival, ";\n",
             "TIMECENSORED = ", mod$input$variable$timecensored, ";\n",
             "CLASSES =  ",mod$input$variable$classes, ";"),
    ANALYSIS = 
      paste0("TYPE = ",mod$input$analysis$type,";\n",
             "DISTRIBUTION = ", mod$input$analysis$distribution, ";\n",
             "PROCESSORS = ", mod$input$analysis$processors, ";\n",
             "STARTS = ", mod$input$analysis$starts, ";\n",
             "STSCALE = ", mod$input$analysis$stscale,";"),
    DEFINE = paste0(mod$input$define,collapse="\n"),
    MODEL =  paste0(mod$input$model, collapse="\n"),
    OUTPUT = paste0(mod$input$output, collapse="\n")
  )
  
  #Dry run to create dataout in desired directory
  mplusModeler(obj, modelout=out,
               dataout = dataout,
               writeData='ifmissing', #This is the only possible option here; set explicitly to suppress a warning message
               hashfilename = hashfilename,
               run=0L
  )
  
  #Once dryrun dataset is made, stringBreak the file.path with custom function
  dataout <- stringBreak(string=dataout,sep="/",accomHash = hashfilename)
  
  return(list(
    obj=obj,
    dataout=dataout
  ))
  
}
