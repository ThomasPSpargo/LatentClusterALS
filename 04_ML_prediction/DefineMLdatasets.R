suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(data.table)
  library(caret)
})
option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file from which datasets will be derived")
  )
opt = parse_args(OptionParser(option_list=option_list))

#Generate directory for output files
path <- file.path(dirname(opt$resultsfile),paste0("MLprediction"))
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import primary dataset and format
data <- tibble(data.table::fread(file=opt$resultsfile))

#Set the reference level as the largest class
tab<- table(data$C)
data$C<- relevel(as.factor(data$C),ref=which(tab==max(tab)))

#Run each each RFE repeating with only people with non-censored disease duration
dfs<- list()

# ### Model with all LCA variables
dfs$data_LCAvars <- data %>%
  dplyr::select("C","CTYDEL","Sex_at_birth","Site_of_Onset","Age_at_onset_years","Time_to_death_or_last_followup_years","Phenotype")

# ### Model with all LCA variables excluding those who are censored 
dfs$data_LCAvars_noncens <- data %>%
  filter(survival_status_bin==1) %>%
  dplyr::select("C","CTYDEL","Sex_at_birth","Site_of_Onset","Age_at_onset_years","Time_to_death_or_last_followup_years","Phenotype")

### Prediction from Pheno vars without considering genetics
dfs$data_pheno <- data %>%
  dplyr::select("C","CTYDEL","Sex_at_birth","Site_of_Onset","Age_at_onset_years","Phenotype")

dfs$data_rarepathways_allPRS <- data %>%
  dplyr::select("C","CTYDEL","Sex_at_birth","Site_of_Onset","Age_at_onset_years","Phenotype","AutoProt","RNAFunc","CytoTransp","c9_status","SOD1","ALScc_SBayesR","FTDcc_SBayesR","SZcc_SBayesR","AZcc_SBayesR","PDcc_SBayesR")

dfs<- lapply(dfs,na.omit)

# ### Create a pheno-data only dataset matched to genetic cohort
dfs$phenofilter_rarepathways_allPRS <- dfs$data_rarepathways_allPRS %>%
  select("C","CTYDEL","Sex_at_birth","Site_of_Onset","Age_at_onset_years","Phenotype")

#### Investigate prediction after filtering down to classes 1 and 2
dfs$class12_pheno <- dfs$data_pheno %>%
  filter(C %in% c(1,2))

dfs$class12_rarepathways_allPRS <- dfs$data_rarepathways_allPRS %>%
  filter(C %in% c(1,2))

dfs$phenofilter_class12_rarepathways_allPRS <- dfs$class12_rarepathways_allPRS %>%
  select("C","CTYDEL","Sex_at_birth","Site_of_Onset","Age_at_onset_years","Phenotype")

#Pivot Phenotype and origin into dummy based on available values in subset.
dfs<- lapply(dfs,function(df){
  
  # One-hot encode (fullrank true for normal dummy coding) relevant columns that may be present within the data
  toDummy <- c("origin","Phenotype")
  toDummy <- toDummy[toDummy %in% names(df)]
  
  if(length(toDummy)>0){
    #Ensure cols are factors and then dummy code
    df[names(df) %in% toDummy] <- df[names(df) %in% toDummy] %>%
      mutate_if(~!is.factor(.),as.factor)
    df <- cbind(df[,!names(df) %in% toDummy],
                predict(dummyVars(reformulate(toDummy),sep="_",df,fullRank = FALSE),df)
    )
  }
  return(df)
})

######
### Write datasets to RDS files
######

for(i in 1:length(dfs)){
  
  #Drop very small factor levels, those with < 20 observations in training data
  checkLevs<- table(dfs[[i]]$C)
  if(any(checkLevs)<20){
    dfs[[i]]<- dfs[[i]] %>%
      filter(!C %in% as.numeric(names(checkLevs)[checkLevs<20]))
  }
  
  #Ensure no empty factor levels and prefix levels with X to ensure no issues parsing in ML predictions
  dfs[[i]]$C <- droplevels(dfs[[i]]$C)
  levels(dfs[[i]]$C) <- paste0("X",levels(dfs[[i]]$C))
  
  #Only write the Dataset file if one isn't already present in the outbound directory
  outfile <- file.path(path,paste0(names(dfs)[i],"_dataset.Rds"))
  if(!file.exists(outfile)){
    saveRDS(dfs[[i]],file=outfile)
  } else{
    message(basename(outfile), " dataset file not written because the file already exists.")
  }
}
