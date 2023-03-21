library(dplyr)
#################################################################################
#Run this script to extract PRS files for neurological traits
options(echo=TRUE)

## Declare path in which the summary file will be returned
outpath <- "/scratch/users/k1802739/GPP_ALS_ProjectMinE/extract_for_LCA"

####
## Read in external files
####

#Read-in stratum 6 IDs to subset ALS prs to only peoople certain not to overlapp with GWAS
s6ID  <- read.delim("/scratch/users/k1802739/GPP_ALS_ProjectMinE/datafiles/PMINE_S6IDs.csv",
                    stringsAsFactors = FALSE,sep = "\t")


#####
## READ IN SBayesR PRS FILES
####

traits <- c("ALScc","FTDposMatch","SZposMatch","AZposMatch","PDposMatch")


root <- "/scratch/users/k1802739/GPP_ALS_ProjectMinE/ALS_ProjectMinE/prs/EUR/"  
list.trait.prs <- list()
for (i in 1:length(traits)) {
  
  #Read sbayesr PRS file
  trait_sbayesr <- read.delim(paste0(root,"sbayesr/",traits[i],"/ALS_ProjectMinE.",traits[i],".EUR.profiles"),sep=" ")
  trait_sbayesr$FID <- NULL #Drop the FID column to avoid issues with merging
  colnames(trait_sbayesr)[1] <- "Illumina_Sample_Barcode" #Change the IID column to match the name of the phenotype ID column
  colnames(trait_sbayesr) <-  gsub("posMatch","cc",colnames(trait_sbayesr)) #Assign name indicating this is a case-control GWAS
  
  #Save in list
  list.trait.prs[[i]] <- trait_sbayesr
  
  names(list.trait.prs)[i] <- gsub("posMatch","cc",traits[i])
}

##### 
# Harmonise into a summary dataframe
####

#Drop any people not in ALS 2022 GWAS Stratum 6 for the ALS polygenic score
list.trait.prs$ALScc<- list.trait.prs$ALScc[list.trait.prs$ALScc$Illumina_Sample_Barcode %in% s6ID$FID2,]

#Compile all secondary traits into one df
prs.all <- list.trait.prs$ALScc %>%
  full_join(list.trait.prs[[2]], by = "Illumina_Sample_Barcode") %>%    
  full_join(list.trait.prs[[3]], by = "Illumina_Sample_Barcode") %>%    
  full_join(list.trait.prs[[4]],  by = "Illumina_Sample_Barcode") %>%
  full_join(list.trait.prs[[5]], by = "Illumina_Sample_Barcode")


##### 
# Write files
####

#Write out file containing all SBayesR PRS
  write.table(prs.all,
              file=file.path(outpath,"SBayesR_PRS_final.csv"),
              sep=",",
              col.names = TRUE,row.names=FALSE)



