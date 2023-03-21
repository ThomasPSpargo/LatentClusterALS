library(dplyr)
library(data.table)

#Specify paths for the input and output files
path <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets"

#Mplus files are specifically directed to a unique directory
mpluspath <- file.path(path,"mplusfiles")
if(!dir.exists(mpluspath)){dir.create(mpluspath)}

# Read in the diagnostic delay standardisation panel
delay_panel<- read.delim(file.path(path,"diagDelay_panel.tsv"),
                         sep="\t", header=TRUE,na.strings = ".")

#Read in PM and STRENGTH datasets
pmAllfields<- fread(file.path(path,"PM_LCA_Allfields.tsv"))
strAllfields<- fread(file.path(path,"STR_LCA_Allfields.tsv"))

#Harmonise STRENGTH column names to correspond with Project MinE column names
strHarmonised <-strAllfields %>%
  rename(Illumina_Sample_Barcode_or_id_orig=id_orig,
         Sex_at_birth=sex,
         Phenotype=diagnosis,
         Site_of_Onset=bulbar,
         Date_of_Onset=date_on,
         Date_of_Diagnosis=date_diag,
         diagnostic_delay_months=diag_delay,
         Date_of_Death=date_death,
         Date_of_survival_check=last_seen,
         Time_from_diagnosis_to_death_or_followup_years=duration,
         survival_status_bin=status,
         Age_at_onset_years=age_onset,
         Diagnostic_delay_years=DELAY,
         Time_to_death_or_last_followup_years=SRV) %>%
  dplyr::select(.,-c(Date_of_Onset,Date_of_Diagnosis,Date_of_Death,Date_of_survival_check))
  
#Harmonise PM ID column names to correspond with (new) STRENGTH ID column
pmHarmonised <- pmAllfields %>%
  rename(Illumina_Sample_Barcode_or_id_orig=Illumina_Sample_Barcode)

#Combine into a single dataset
fulldsetAllfields<- full_join(pmHarmonised,strHarmonised)
#nrow(pmAllfields)+nrow(strHarmonised) #14352
#nrow(fulldsetAllfields) #14352

######
### Generate per-country Standardised diagnostic delay column
######

#Define function for standardising diagnostic delay
scaleToPanel <- function(x,ID,panel){
  rowID<- panel$ID==ID        #define the panel row to scale against
  x <- x-panel$mean[rowID]    #center on mean
  x <- x/panel$sd[rowID]      #Scale to reference variance 
  
  x
}

#Scale country wise diagnostic delay according to panel
countries<- unique(fulldsetAllfields$origin)
fulldsetAllfields$CTYDEL <- NA
for(i in 1:length(countries)){
  origin.ind<- fulldsetAllfields$origin==countries[i]
  fulldsetAllfields$CTYDEL[origin.ind] <- scaleToPanel(fulldsetAllfields$Diagnostic_delay_years[origin.ind],ID = countries[i],panel=delay_panel)
}


### Write to file the final dataset before converting to Mplus format
write.table(fulldsetAllfields,file=file.path(path,"JNT_LCA_Allfields.tsv"),row.names=FALSE,col.names = TRUE,sep="\t",quote = FALSE)

# #### Update the PRS
# 
# newPRS<- read.csv("/Users/tom/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/datasets/SBayesR_PRS_final.csv")
# allFields<- tibble(fread(file.path(path,"JNT_LCA_Allfields.tsv")))
# 
# #Assign matching name for ID column
# colnames(newPRS)[colnames(newPRS)=="Illumina_Sample_Barcode"] <- "Illumina_Sample_Barcode_or_id_orig"
# 
# #Check that matches are all to PM
# # length(which(newPRS$Illumina_Sample_Barcode_or_id_orig %in% allFields$Illumina_Sample_Barcode_or_id_orig[grepl("^999", allFields$new_num_id)]))
# # length(newPRS$Illumina_Sample_Barcode_or_id_orig %in% allFields$Illumina_Sample_Barcode_or_id_orig)
# 
# # #Drop old PRS fields
# ncol(allFields) #108
# allFields<- allFields[,!grepl("SBayesR",colnames(allFields))]
# ncol(allFields) #103
# 
# #Reinsert the DF and write as above
# fulldsetAllfields <- left_join(allFields,newPRS)

#Repeat for controls
# outpath <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets"
# 
# newPRS<- read.csv("/Users/tom/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/datasets/SBayesR_PRS_final.csv")
# CTRLFields<- tibble(fread("/Users/tom/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets/PM_CONTROLS_data.tsv"))
# 
# # #Drop old PRS fields
# ncol(CTRLFields) #42
# CTRLFields<- CTRLFields[,!grepl("SBayesR",colnames(CTRLFields))]
# ncol(CTRLFields) #37
# 
# #Reinsert the DF and write as above
# PM_control <- left_join(CTRLFields,newPRS)
# 
# write.table(PM_control,file=file.path(outpath,"PM_CONTROLS_data.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)



#######
### Prepare input files for running mplus
#######

##### Generate Main LCA sample file - filtering down to the required columns
Mplus_full <- tibble(fulldsetAllfields)[c("new_num_id",
                            "Phenotype",
                            "Sex_at_birth",
                            "Age_at_onset_years",
                            "CTYDEL",
                            "Time_to_death_or_last_followup_years",
                            "survival_status_bin",
                            "Site_of_Onset"
)]

#Rename selected colnames into MPLUS friendly form variables
colnames(Mplus_full) <- c("ID", "PHENO", "SEX", "AGEONS", "CTYDEL", "SRV", "SRVBIN", "SITBIN")


#### Generate the discovery and validation sample data
Mplus_pm<- Mplus_full %>%
  filter(grepl("^999", ID))

####
Mplus_str<- Mplus_full %>%
  filter(grepl("^555", ID))


# mplus cannot handle columns with names. Therefore write to files custom wrapper function around write and write.table functionality
#Dset is the dataset to write
#pathprefix is the file path and the filename for the dataset omitting any extension - it will be written to a .csv file and colnames .txt file 
writeMplus <- function(dset,pathprefix){
  
  #Write the dataset
  write.table(dset,file=paste0(pathprefix,".csv"),row.names=FALSE,na = ".",col.names = FALSE,sep=",",quote = FALSE)
  #Write the column names
  write(colnames(dset),file=paste0(pathprefix,"_colnames.txt"))
}

#Write the full joint, discovery, and validation datasets
writeMplus(Mplus_full,file.path(mpluspath,"JNT_FULL_ANY"))
writeMplus(Mplus_pm,file.path(mpluspath,"PM_FULL_ANY"))
writeMplus(Mplus_str,file.path(mpluspath,"STR_NOPM_ANY"))


### Write the datasets with missingness in Diag delay and disease duration excluded
Mplus_full %>%
  filter(!is.na(CTYDEL) & !is.na(SRV)) %>%
writeMplus(.,file.path(mpluspath,"JNT_FULL_DDOMIT"))

Mplus_pm %>%
  filter(!is.na(CTYDEL) & !is.na(SRV)) %>%
writeMplus(.,file.path(mpluspath,"PM_FULL_DDOMIT"))

Mplus_str %>%
  filter(!is.na(CTYDEL) & !is.na(SRV)) %>%
writeMplus(.,file.path(mpluspath,"STR_NOPM_DDOMIT"))



### Write the files excluding people with ANY missingness
Mplus_full %>%
  na.omit() %>%
  writeMplus(.,file.path(mpluspath,"JNT_FULL_NAOMIT"))

Mplus_pm %>%
  na.omit() %>%
  writeMplus(.,file.path(mpluspath,"PM_FULL_NAOMIT"))

Mplus_str %>%
  na.omit() %>%
  writeMplus(.,file.path(mpluspath,"STR_NOPM_NAOMIT"))

