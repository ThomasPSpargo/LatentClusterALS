
suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(survival)
  library(optparse)
})


option_list = list(
  make_option("--filetoparse", action="store", default=NA, type='character',
              help="Specify the source file for which standardised columns will be generated"),
  make_option("--normalisationPanel", action="store", default=NULL, type='character',
              help="Optional file containing mean and SD for the variables: 'Age_at_onset_years' 'Diagnostic_delay_years' 'Time_to_death_or_last_followup_years'"),
  make_option("--writePanel", action="store", default=FALSE, type='logical',
              help="Logical; specify whether to write a file containing normalisation details for 'Age_at_onset_years' 'Diagnostic_delay_years' 'Time_to_death_or_last_followup_years' in the base data. is ignored if normalisationPanel is specified as input."),
  make_option("--writeSrvSubset", action="store", default=FALSE, type='logical',
              help="Write files into two separate subdirectories, differentiating between the complete dataset and a dataset which omits people with censored disease duration."),
  make_option("--reorderClasses", action="store", default=FALSE, type='logical',
              help="If TRUE, rearrange the class levels such that class 1 is the class with the biggest N and the others are arranged descending")
  
)

opt = parse_args(OptionParser(option_list=option_list))

#Take the name of the dataset specifically
name <- gsub(".tsv","",basename(opt$filetoparse))

path <- dirname(opt$filetoparse)
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import input datafiles
data <- tibble(fread(opt$filetoparse))

if(opt$reorderClasses){
  cat("Classes reordered, pre-format:\n")
  print(table(data$C))
  
  data$C<- as.factor(data$C) #Store as factor
  class_size<- table(data$C) %>% #Determine class order from largest to smallest N and store index position
    sort() %>%
    rev() %>%
    names() %>%
    as.numeric()
  
  levels(data$C)[class_size] <- 1:length(levels(data$C)) #Relevel classes according to class order
  names(data)[which(names(data) %in% paste0("CPROB",class_size))[class_size]] <- paste0("CPROB",1:length(levels(data$C))) #Reorder CPROBS equivalently

  cat("Classes reordered, post-format:\n")
  print(table(data$C))
}

if(!is.null(opt$normalisationPanel)){
  #Read-in the panel
  panel <- fread(opt$normalisationPanel)
  
  #Define function to scale a continuous variable according to mean and SD values indicated in a reference panel
  #x is a numeric vector
  #ID is an identifier for the row against which to scale
  #panel defines the reference panel, with at least the columns: "ID",  "mean",  "sd"
  scaleToPanel <- function(x,ID,panel){
    rowID<- panel$ID==ID        #define the panel row to scale against
    x <- x-panel$mean[rowID]    #center on mean
    x <- x/panel$sd[rowID]      #Scale to reference variance 
    return(x)
  }
  
} else {
  if(opt$writePanel){
    panelpath<-file.path(path,"standardisationPanel.tsv")
    write("ID\tmean\tsd",file=panelpath) #Empty tsv file containing headers for standardisation variables
  }
}

######
### Rescale and centre continuous variables to have Mean=0 and SD=1
### Optionally base these on an existing normalisation panel
### Optionally write out the standardisation scores
######
if("Age_at_onset_years" %in% colnames(data)) {
  if(!is.null(opt$normalisationPanel)){
    data$AGEONS_nml <- scaleToPanel(data$Age_at_onset_years,"Age_at_onset_years",panel)
  } else {
    data$AGEONS_nml <- as.numeric(scale(data$Age_at_onset_years))
    if(opt$writePanel){
      write(paste("Age_at_onset_years",
                  mean(data$Age_at_onset_years,na.rm=TRUE),
                  sd(data$Age_at_onset_years,na.rm=TRUE),
                  sep="\t"),file=panelpath,append=TRUE)
    }
  }
}
if("Diagnostic_delay_years" %in% colnames(data)) {
  if(!is.null(opt$normalisationPanel)){
    data$DELAY_nml <- scaleToPanel(data$Diagnostic_delay_years,"Diagnostic_delay_years",panel)
  } else {
    data$DELAY_nml <- as.numeric(scale(data$Diagnostic_delay_years))
    if(opt$writePanel){
      write(paste("Diagnostic_delay_years",
                  mean(data$Diagnostic_delay_years,na.rm=TRUE),
                  sd(data$Diagnostic_delay_years,na.rm=TRUE),
                  sep="\t"),file=panelpath,append=TRUE)
    }
  }
}
if("Time_to_death_or_last_followup_years" %in% colnames(data)) {
  if(!is.null(opt$normalisationPanel)){
    data$SRV_nml <- scaleToPanel(data$Time_to_death_or_last_followup_years,"Time_to_death_or_last_followup_years",panel)
  } else {
    data$SRV_nml <- as.numeric(scale(data$Time_to_death_or_last_followup_years))
    if(opt$writePanel){
      write(paste("Time_to_death_or_last_followup_years",
                  mean(data$Time_to_death_or_last_followup_years,na.rm=TRUE),
                  sd(data$Time_to_death_or_last_followup_years,na.rm=TRUE),
                  sep="\t"),file=panelpath,append=TRUE)
    }
  }
}
######
### Optionally write files for the complete sample and after removing people with censored survival. If this is done files will be output to a subdirectory
######
if(opt$writeSrvSubset){
  #Write files for subset dataset
  data_subset <- data[which(data$survival_status_bin==1),]
  noncensored_subpath<- file.path(dirname(opt$filetoparse),"noncensored")
  if(!dir.exists(noncensored_subpath)){dir.create(noncensored_subpath,recursive = TRUE)}

  fwrite(data_subset,file=file.path(noncensored_subpath,basename(opt$filetoparse)),sep="\t")
  
  #write files for 'complete' dataset
  fulldata_subpath<- file.path(dirname(opt$filetoparse),"completesample")
  if(!dir.exists(fulldata_subpath)){dir.create(fulldata_subpath,recursive = TRUE)}
  fwrite(data,file=file.path(fulldata_subpath,basename(opt$filetoparse)),sep="\t")
  
} else {
  fwrite(data,file=opt$filetoparse,sep="\t")
}

