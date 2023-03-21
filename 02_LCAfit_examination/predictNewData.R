suppressPackageStartupMessages({
  library(MplusAutomation)
  library(optparse)
  library(data.table)
})

option_list = list(
  make_option("--SOURCEFILE", action="store", default=NA, type='character',
              help="Optionally specify the source file from which a model will be obtained. This file must have the 'SVALUES' option of Mplus set, in order to determine parameter configuration. If option is left blank, the parameters from the accepted 5-class joint dataset LCA model associated with this repository will be used as a basis."),
  make_option("--DATASET", action="store", default=NA, type='character',
              help="Specify path to dataset to which model will be applied, expects the format: <file>.csv and a corresponding file with column names: <file>_colnames.txt"),
  make_option("--USEOBSERVATIONS", action="store", default=NA, type='character',
              help="Indicate criteria for data filtering (see USEOBSERVATIONS in mplus"),
  make_option("--DEFINE", action="store", default=NA, type='character',
              help="Indicate variable mutate options (see DEFINE in Mplus)"),
  make_option("--out", action="store", default=NA, type='character',
              help="Indicate the directory and prefix for output directory"),
  make_option("--PreMplusFields", action="store", default=NA, type='character',
              help="Specify a dataset to merge with the MPLUS output, taking minimal columns from MPlus (Cprobs and class assignment)"),
  make_option("--PreMplusIDname", action="store", default="ID", type='character',
              help="Specify the column name in the pre mplus dataset upon which to match, defaults to 'ID'"),
  make_option("--MplusKeepcolumns", action="store", default="*", type='character',
              help="Pattern passed to grepl, indicates which columns from the Mplus results should be retained. Defaults to '*', which returns all columns. Use '|' to indicate OR for multiple patterns"),
  make_option("--helperFunDir", action="store", default=".", type='character',
              help="Directory containing custom helper functions leveraged in the script"),
  make_option("--forceNoDummies", action="store", default="FALSE", type='logical',
              help="Logical, set to TRUE to disable the creation of dummy records in the new dataset to 'fill-in' missing categorical levels when passing to mplus. These records are only created if certain categorical variable levels are entirely missing from the new dataset, and are removed when the data is returned to the user. When levels are missing Mplus may not allow the model to be applied to the new data. This option is provided in case users would not want dummys in their data under any circumstances but it may cause the script to fail.")
)

opt = parse_args(OptionParser(option_list=option_list))

#Read-in all helper functions
invisible(lapply(list.files(opt$helperFunDir,full.names = TRUE,pattern=".R"),source))

######
### Prepare inputs to pass to MPLUS according to options set in useobservations and define
######
#If any USEOBSERVATIONS or DEFINE options are specified, split the comma-delimited string and then collapse by semicolon and newline
appendout <- "" #Additionally append relevant string to the end of the directory path to declare options
if(!is.na(opt$USEOBSERVATIONS)){
  useobs <- paste(strsplit(opt$USEOBSERVATIONS,split=",")[[1]],collapse=";\n")
  
  appendout <- paste(appendout,"withdrop",sep = "_")
} 
if(!is.na(opt$DEFINE)){
  defineobs <- paste(strsplit(opt$DEFINE,split=",")[[1]],collapse=";\n")
  appendout <- paste(appendout,"withmiss",sep = "_")
}

#If there is no modification, return as-is
if(appendout=="") appendout="_as_is"

#Import new dataset and assign column names (if none detected)
newdata <- fread(file=opt$DATASET, na.strings = ".")

#If all column names are default, assign from matched column-name file
if(all(grepl("V",names(newdata)))){
  newdata_colnames <- gsub(".csv","_colnames.txt", opt$DATASET)
  names(newdata)<- scan(newdata_colnames,what="character")
}

#If there is a user-supplied model, determine the input for the new data from this
if(!is.na(opt$SOURCEFILE)){
  
  #Extract fixed model outputs from the supplied SOURCEFILE
  path<- opt$SOURCEFILE
  
  #Read in the original fit, but skip the savedata and bparameters options
  mod <- readModels(path,
                    what=c("input", "warn_err", "data_summary", "sampstat", "covariance_coverage", "summaries",
                           "parameters", "class_counts", "indirect", "mod_indices", "residuals",
                           "tech1", "tech3", "tech4", "tech7", "tech8", "tech9", "tech10", "tech12",
                           "fac_score_stats", "lcCondMeans", "gh5", "output"))
  
  #Identify the start and endpoint for the SVALUES output
  commandstarts<- grep("MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES", mod$output)+1
  commandends<- grep("TECHNICAL 1 OUTPUT", mod$output)-1
  
  #Extract the results and fix param values
  fixparams <- paste0(mod$output[commandstarts:commandends],collapse="\n")
  fixparams <- gsub("\\*","@",fixparams)
  
  #Extract posterior probabilities of each class
  nclasses<- mod$input$variable$classes
  
  #Store in list
  fixed <- list(fixparams=fixparams, nclasses=nclasses, mod = mod)
  
} else {
  #Use the pre-defined model, extracted from the final LCA fit applied to the joint dataset. Object structured to mimic MplusAutomation object structure
  fixed <- list(
    fixparams="\n     %OVERALL%\n\n     [ c#1@-1.67289 ];\n     [ c#2@-4.24706 ];\n     [ c#3@-0.22669 ];\n     [ c#4@-3.14519 ];\n\n     %C#1%\n\n     [ ageons@58.77206 ];\n     [ ctydel@1.36767 ];\n     [ srv@-4.33335 ];\n     [ pheno#1@2.41551 ];\n     [ pheno#2@0.61141 ];\n\n     [ sex$1@0.39001 ];\n     [ sitbin$1@1.81760 ];\n\n     ageons@134.19644 (3);\n     ctydel@0.17535 (4);\n\n     %C#2%\n\n     [ ageons@40.05059 ];\n     [ ctydel@8.08113 ];\n     [ srv@-5.96121 ];\n     [ pheno#1@2.01613 ];\n     [ pheno#2@1.03181 ];\n\n     [ sex$1@0.65893 ];\n     [ sitbin$1@2.09112 ];\n\n     ageons@134.19644 (3);\n     ctydel@0.17535 (4);\n\n     %C#3%\n\n     [ ageons@57.46387 ];\n     [ ctydel@-0.11572 ];\n     [ srv@-3.64307 ];\n     [ pheno#1@2.49259 ];\n     [ pheno#2@-0.60018 ];\n\n     [ sex$1@0.53478 ];\n     [ sitbin$1@1.40931 ];\n\n     ageons@134.19644 (3);\n     ctydel@0.17535 (4);\n\n     %C#4%\n\n     [ ageons@55.50137 ];\n     [ ctydel@3.62096 ];\n     [ srv@-5.07112 ];\n     [ pheno#1@1.93772 ];\n     [ pheno#2@1.17556 ];\n\n     [ sex$1@0.49973 ];\n     [ sitbin$1@1.75552 ];\n\n     ageons@134.19644 (3);\n     ctydel@0.17535 (4);\n\n     %C#5%\n\n     [ ageons@64.65670 ];\n     [ ctydel@-0.38158 ];\n     [ srv@0 ];\n     [ pheno#1@3.73459 ];\n     [ pheno#2@-2.01095 ];\n\n     [ sex$1@0.25042 ];\n     [ sitbin$1@0.34671 ];\n\n     ageons@134.19644 (3);\n     ctydel@0.17535 (4);\n\n\n",
    nclasses="c (5)",
    mod = list(
      input=list(
        title="JNT_FULL_ANY_LCA_c5_savedata",
        variable=list(
          usevariables="ID PHENO SEX AGEONS CTYDEL SITBIN SRV SRVBIN",
          useobservations=NULL,
          idvariable="ID",
          nominal="PHENO",
          categorical="SEX SITBIN",
          survival="SRV",
          timecensored="SRVBIN (1 = NOT 0 = RIGHT)",
          classes="c (5)"),
        analysis=list(
          type="MIXTURE",
          distribution="",
          processors="8",
          starts="0",
          stscale="15"),
        define=NULL,
        model=NULL,
        output=c("","","TECH1;","ENTROPY;","SVALUES;")
      )
    )
  )
}

#Function to apply fixed model to secondary dataset 
#Fixed - extraction of features from training dataset for application to new dataset
#newdata - new dataset containing at least the same features as the usevariables of the first dataset
#opt - settings from optparse when running the script
predictNew <- function(fixed,newdata,opt){
  
  #Prepare the output directory
  out_dir <- paste0(opt$out,"_predict",appendout)
  if(!dir.exists(out_dir)){dir.create(out_dir,recursive = TRUE)}
  
  #Identify class number in training dataset
  m<- regexpr("[0-9]+", fixed$nclasses)
  numclass<- regmatches(fixed$nclasses,m)
  
  #Make the file output location - name according to numclass
  out<- file.path(out_dir,"LCAnewdata.inp")
  dataout <- file.path(out_dir,"LCAnewdata.dat")
  
  #When there are categorical variables with levels missing in the dataset, mplus can return an error message.
  #Unless --forceNoDummies is TRUE, fill in the blank categorical levels with dummy records
  if(!opt$forceNoDummies && is.na(opt$SOURCEFILE)){
    #Determine if there are any empty levels in the data
    phenomissing<- which(!(c(1:3) %in% sort(unique(newdata$PHENO))))
    sexmissing <- which(!(c(1:2) %in% sort(unique(newdata$SEX))))
    sitbinmissing <- which(!(c(1:2) %in% sort(unique(newdata$SITBIN))))
    if(sum(phenomissing,sexmissing,sitbinmissing)>0){
      dummyID1<- dummyID<- 8888880
      
      #Phenomissing indicates missing categorical variable levels
      #newdata is the secondary dataset, and new rows are added for the missing people
      if(length(phenomissing)>0){
        for(i in 1:length(phenomissing)){
          newdata <- rbind(newdata, list(ID=dummyID,PHENO=phenomissing[i],SRVBIN=1),fill=TRUE)
          dummyID <- dummyID+1
        }
      }
      if(length(sexmissing)>0){
        for(i in 1:length(sexmissing)){
          newdata <- rbind(newdata, list(ID=dummyID,SEX=sexmissing[i],SRVBIN=1),fill=TRUE)
          dummyID <- dummyID+1
        }
      }
      if(length(sitbinmissing)>0){
        for(i in 1:length(sexmissing)){
          newdata <- rbind(newdata, list(ID=dummyID,SITBIN=sitbinmissing[i],SRVBIN=1),fill=TRUE)
          dummyID <- dummyID+1
        }
      }
      warning("Certain levels were missing from categorical variables in the new dataset. This can cause issues with Mplus. Accordingly, dummy records were temporarily appended to the dataset, using IDs:",paste0(dummyID1:dummyID, collapse=", "),". They have been again removed from the final dataset returned to the user.\n\n")
      finalDummies<- dummyID1:dummyID 
    }
  } else {
    finalDummies <- NULL
  }
  
  #Convert the model into a mplusAutomation mplusObject with custom function and write the datafile to the output directory
  #Filename is not hashed to avoid issues with dropping people from the dataset
  prep <- prepObject(fixed$mod,out,dataout,newData = newdata,hashfilename=FALSE)
  obj <- prep$obj
  dataout <- prep$dataout
  
  #Apply a title 
  obj$TITLE <- "LCAfit_new_data"
  
  #Set to 0 starts
  obj$ANALYSIS <- gsub(paste0("STARTS = +[0-9]+( +[0-9]+;|;)"), 
                       paste0("STARTS = 0;"),
                       obj$ANALYSIS)
  
  #Assign savedata params, defining the name for the output file
  obj$SAVEDATA <- paste0("FILE IS ", obj$TITLE,".dat;\n",
                         "FORMAT IS FREE;\n",
                         "SAVE = CPROB;\n") #Set option to save class probability data  
  
  #If useobservations is set, append this call
  if(!is.na(opt$USEOBSERVATIONS)){
    if(grepl("USEOBSERVATIONS =",obj$VARIABLE)) warning("USE OBSERVATIONS ALREADY DEFINED PLEASE CHECK THAT ERRONEOUS SCRIPTS HAVE NOT BEEN INTRODUCED")
    
    obj$VARIABLE <- paste0(obj$VARIABLE,"\n","USEOBSERVATIONS = ", useobs,";")
  }
  
  if(!is.na(opt$DEFINE)) obj$DEFINE <- paste0(defineobs,";")
  
  ### If the model parameters are absent in the new data, check, adjust accordingly, and return warning
  vars<- strsplit(fixed$mod$input$variable$usevariables,split=" ")[[1]]
  vars <- vars[vars!=""]
  mis_col<- vars[which(!(vars %in% colnames(newdata)))]
  
  if(length(mis_col)>0){
    warning("Mplus expects the following input usevariables: ", paste(vars,collapse = " "),".\n",
            "The new dataset provided is missing the following: ",paste(mis_col,collapse= " "),".\n",
            "The model command has been adjusted to exclude the missing columns; please check that this is acceptable in your use-case."
    )
    
    params <- strsplit(fixed$fixparams,split="\n")[[1]]
    #String split and drop params referring to missing variables, drop also from the variable command
    for(i in 1:length(mis_col)){
      params <-params[!grepl(paste0("(?i)",tolower(mis_col[i]),"(#|\\;|\\@)"),params)]
      obj$VARIABLE <- gsub(mis_col[i],"",obj$VARIABLE)
    }
    #Recombine the params and add to model command
    obj$MODEL <- paste0(params,collapse="\n")
    
  } else {
    obj$MODEL <- fixed$fixparams  #If all cols are present, add to fixed model parameters to model command directly
  }

  #Real run, with line-break for dataout to avoid issue with mplus reading the file
  mod2 <- mplusModeler(obj, modelout=out,
                       dataout = dataout,
                       writeData='never', #This is the only allowed option; set to suppress a warning message
                       hashfilename = FALSE, #FALSE in order to avoid 90 character file limits
                       run=1L)
  
  #Extract the mplus-saved data
  data <- mod2$results$savedata
  
  #If dummy records were created, remove them before returning the dataset to the user
  if(!is.null(finalDummies)) data <- data[!data$ID %in% finalDummies,]

  #If using the PRESET mplus file, reorder classes to correspond with class order from manuscript
  if(is.na(opt$SOURCEFILE)){
    message("NOTE: The numerical order of classes assigned by Mplus has been reassigned to match with the published class order.\n\nClass orders are as follows [mplus-assigned number -> published class number]:\n1 to 3\n2 to 5\n3 to 2\n4 to 4\n5 to 1\n\n")
    
    #Reorder classes and descriptive CPROB columns
    publishOrder<- c(3,5,2,4,1)
    data$C <- factor(data$C,levels=1:5,labels=publishOrder)
    colnames(data)[grepl("CPROB",colnames(data))] <- paste0("CPROB",publishOrder)
  }
  
  
  
  if(!is.na(opt$PreMplusFields)){ #If a dataset has been provided for merging additional cols, prepare and combine the Mplus and extra data
    fullDF<- fread(opt$PreMplusFields)
    
    colnames(data)[colnames(data)=="ID"] <- opt$PreMplusIDname
    keepcols<- c(opt$PreMplusIDname,grep(opt$MplusKeepcolumns,colnames(data),value=TRUE))
    
    message("Joining datasets and Keeping the mplus columns: ",paste(keepcols,collapse= ", "))
    data <- left_join(data[,keepcols],fullDF,by=opt$PreMplusIDname) #Keep only common ID columns
  }
  
  #Write final tsv file
  finalfile <- file.path(out_dir,paste0(obj$TITLE,".tsv"))
  fwrite(data,file=finalfile,sep="\t")
  
  #Remove intermediate dataset files produced by mplus
  system(paste0("rm ",out_dir,"/*.dat"))
  
  message(obj$TITLE," file written")
  
  return(finalfile)
}

outfile<- predictNew(fixed,newdata,opt)

#Run function
cat(outfile)
