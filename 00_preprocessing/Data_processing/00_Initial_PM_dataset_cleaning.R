#rm(list=ls())

library(lubridate) #For days to years conversion
library(dplyr) #%>%
library(readxl) #For reading genetic dataset(s)
library(tidyr) #For pivoting

#Declare path for all datasets
path <- "/Users/tom/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/datasets"


#New Pheno data script
NP<- readRDS(file = file.path(path,"pheno_input.2020-10-29_WGS_onset_KLC.rds"))

#Read old pheno data (includes QC and principal components)
PM_phen <- read.delim(file.path(path,"recodePM.pheno.2_1_correct.txt"))

#Harmonise current and legacy data
# # Check that NP data includes all samples from PM_Phen
# identical(
#   length(PM_phen$Illumina_Sample_Barcode),
#   length(which(PM_phen$Illumina_Sample_Barcode %in% NP$Illumina_Sample_Barcode))
# )
# #[1] TRUE

#Take relevant QC columns from old data
PM_tojoin <- PM_phen %>%
  dplyr::select(Illumina_Sample_Barcode,failed_QC, exclude_reason,relatedness)

#Left-join these columns to the NP dataset
NP <- NP %>%
  left_join(PM_tojoin,by = colnames(NP)[1])       #Binary vars

#See how many LP barcodes lack a QC pass value
# noqcdata<-which(is.na(PMFULLdat$failed_QC))
# length(grep("LP",PMFULLdat$Illumina_Sample_Barcode[noqcdata])) #519

#write.csv(NP,file=paste0("New_PM_Pheno",Sys.Date(),".csv"))

#Remove the old data from environment
rm(list=ls(pattern="PM"))


#Set missing some impossible values:
NP$Age_at_ventilation_more_23h_per_day_days[which(NP$Age_at_onset_days>NP$Age_at_ventilation_more_23h_per_day_days)] <- NA
NP$Date_of_ventilation_more_23h_per_day[which(NP$Date_of_Onset>NP$Date_of_ventilation_more_23h_per_day)] <- NA
NP$Age_at_diagnosis_days[which(NP$Age_at_onset_days>NP$Age_at_diagnosis_days)] <- NA
NP$Date_of_Diagnosis[which(NP$Date_of_Onset>NP$Date_of_Diagnosis)] <- NA

#None
#NP$Age_at_death_days[which(NP$Age_at_onset_days>NP$Age_at_death_days)] <- NA


#--------------------------------------------------------------#
#           Convert continuous variables into years            #
#--------------------------------------------------------------#
##Computing: "Age_at_onset_years"
#Try based on Age_at_onset_days
#If nothing, then try based on interval(NP$Date_of_Birth,NP$Date_of_Onset)

#Attempt to identify value from Age_at_onset_days:
#Round down to nearest whole day
d.floor<- floor(NP$Age_at_onset_days)
#Calculate number of partial days
remainder <- NP$Age_at_onset_days-d.floor
#Convert remainder into hours  
rem.inHours<- remainder*24

#Calculate length of time using lubridate package, to the nearest hour
lub<- days(d.floor)+hours(round(rem.inHours))


#If doesn't work, assign values based on dates of birth and onset
#Index people who have not been assigned a value
empty<- which(is.na(lub))

#Create interval information
#https://stackoverflow.com/questions/63523659/convert-days-to-years-in-r
#https://stackoverflow.com/questions/15569333/get-date-difference-in-years-floating-point #2nd link for time_length function
lub[empty]<-lubridate::as.period(lubridate::interval(NP$Date_of_Birth[empty],NP$Date_of_Onset[empty]))
#length(which(!is.na(lub[empty]))) #0 records updated in this way

#Convert everything into time-length data in year units
NP$Age_at_onset_years <- lubridate::time_length(lub,unit="year")


#------------#

##Computing: "Diagnostic_delay_years"
#Try based on Diagnostic_delay_days
#If nothing, then try Age_at_diagnosis_days-Age_at_onset_days
#If nothing, then try based on interval(NP$Date_of_Diagnosis-NP$Date_of_Onset)

#Round down to nearest whole day
d.floor<- floor(NP$Diagnostic_delay_days)
#Calculate number of partial days
remainder <- NP$Diagnostic_delay_days-d.floor
#Convert remainder into hours  
rem.inHours<- remainder*24

#Calculate length of time using lubridate package, to the nearest hour
lub<- days(d.floor)+hours(round(rem.inHours))


#If doesn't work, assign values based on age of onset and diagnosis
#Index people who have not been assigned a value
empty<- which(is.na(lub))
#length(empty) #9558

#Calculate diagnostic delay
interval <- NP$Age_at_diagnosis_days[empty]-NP$Age_at_onset_days[empty]
#Round down to nearest whole day
d.floor<- floor(interval)
#Calculate number of partial days
remainder <- interval-d.floor
#Convert remainder into hours  
rem.inHours<- remainder*24

#Calculate length of time using lubridate package, to the nearest hour
lub[empty]<- days(d.floor)+hours(round(rem.inHours))
#length(which(!is.na(lub[empty]))) #40 records updated in this way


#If doesn't work, assign values based on date of onset and diagnosis
#Index people who have not been assigned a value
empty<- which(is.na(lub))
#length(empty) #9518

#Create interval information
lub[empty]<-lubridate::as.period(lubridate::interval(NP$Date_of_Onset[empty],NP$Date_of_Diagnosis[empty]))
#length(which(!is.na(lub[empty]))) #0 records updated in this way

#Convert to time length in unit years
NP$Diagnostic_delay_years <- lubridate::time_length(lub,unit="year")


#------------#    

##Computing: "Time_to_death_or_last_followup_years"
#Determine time to death variable:
#First by Age_at_death_days-Age_at_onset_days
#Then if not by interval(NP$Date_of_Onset,NP$Date_of_Death)

#Calculate time to death
interval <- NP$Age_at_death_days-NP$Age_at_onset_days
#Round down to nearest whole day
d.floor<- floor(interval)
#Calculate number of partial days
remainder <- interval-d.floor
#Convert remainder into hours  
rem.inHours<- remainder*24

#Calculate length of time using lubridate package, to the nearest hour
lub<- days(d.floor)+hours(round(rem.inHours))
#length(which(!is.na(lub))) #8015 records updated in this way


#If doesn't work, assign values based on date of onset and death
#Index people who have not been assigned a value
empty<- which(is.na(lub))
#length(empty) #10140

#Create interval information
lub[empty]<-lubridate::as.period(lubridate::interval(NP$Date_of_Onset[empty],NP$Date_of_Death[empty]))
#length(which(!is.na(lub[empty]))) #1 record updated in this way


#Convert to time length in unit years
NP$Time_to_death <- lubridate::time_length(lub,unit="year")

####
#Related to the survival check:
#Attempt first based on interval(NP$Date_of_Onset,NP$Date_of_survival_check)
#If nothing, then try based on Age_at_survival_check_days-Age_at_onset_days
#Date_of_Onset
#Date_of_Death
#Date_of_survival_check
#"Date_of_ventilation_more_23h_per_day"
####
#Determine time to survival check
#First by Age_at_survival_check_days-Age_at_onset_days
#Then if not by interval(NP$Date_of_Onset,NP$Date_of_survival_check)

#Calculate time to death
interval <- NP$Age_at_survival_check_days-NP$Age_at_onset_days
#Round down to nearest whole day
d.floor<- floor(interval)
#Calculate number of partial days
remainder <- interval-d.floor
#Convert remainder into hours  
rem.inHours<- remainder*24

#Calculate length of time using lubridate package, to the nearest hour
lub<- days(d.floor)+hours(round(rem.inHours))
#length(which(!is.na(lub))) #8015 records updated in this way


#If doesn't work, assign values based on date of onset and surv check
#Index people who have not been assigned a value
empty<- which(is.na(lub))
#length(empty) #7775

#Create interval information
lub[empty]<-lubridate::as.period(lubridate::interval(NP$Date_of_Onset[empty],NP$Date_of_survival_check[empty]))
#length(which(!is.na(lub[empty]))) #4 records updated in this way

#Convert to time length in unit years
NP$Time_to_survival_check <- lubridate::time_length(lub,unit="year")



#Determine time to ventilation
#First by Age_at_ventilation_more_23h_per_day_days-Age_at_onset_days
#Then if not by interval(NP$Date_of_Onset,NP$Date_of_ventilation_more_23h_per_day)

#Note that only a small number of people have this variable
#length(which(is.na(NP$Date_of_ventilation_more_23h_per_day)))  #17997
#length(which(!is.na(NP$Date_of_ventilation_more_23h_per_day))) #158

#Calculate time to death
interval <- NP$Age_at_ventilation_more_23h_per_day_days-NP$Age_at_onset_days
#Round down to nearest whole day
d.floor<- floor(interval)
#Calculate number of partial days
remainder <- interval-d.floor
#Convert remainder into hours  
rem.inHours<- remainder*24

#Calculate length of time using lubridate package, to the nearest hour
lub<- days(d.floor)+hours(round(rem.inHours))
#length(which(!is.na(lub))) #168 records updated in this way


#If doesn't work, assign values based on date of onset and ventilation
#Index people who have not been assigned a value
empty<- which(is.na(lub))
#length(empty) #17987

#Create interval information
lub[empty]<-lubridate::as.period(lubridate::interval(NP$Date_of_Onset[empty],NP$Date_of_ventilation_more_23h_per_day[empty]))
#length(which(!is.na(lub[empty]))) #0 records updated in this way



#Convert to time length in unit years
NP$Time_to_ventilation <- lubridate::time_length(lub,unit="year")



#### Compare survival check and death data
#which(NP$Time_to_survival_check > NP$Time_to_death)                   #integer(0) - no survival check is recorded after death
#which(NP$Time_to_survival_check < NP$Time_to_death)                   #integer(0) - no deaths are recorded later than survival check
#which( is.na(NP$Time_to_survival_check) & !is.na(NP$Time_to_death))   #integer(0) - No people are NA on survival check and have TTD data
#which(!is.na(NP$Time_to_survival_check) &  is.na(NP$Time_to_death))   #             Many people are NA on TTD and have Survival check data (expected, as these will be censored)

#Identify people with ventilation >23h/day, and survival date passes their ventilation date
VentAsDeath<-which(NP$Time_to_survival_check > NP$Time_to_ventilation)  #Many people have survived longer than TTV
#which(NP$Time_to_survival_check < NP$Time_to_ventilation)             #3 People survival check shorter than time to ventilation

#ventBeforeDeath<- which(NP$Time_to_death>NP$Time_to_ventilation) #(expected)
#ventAfterDeath<- which(NP$Time_to_death<NP$Time_to_ventilation) #3 people ventilated after death; considered deceased as in time-to-death

#Compute final survival variable
#Since time_to_death is captured properly in the survival check variable, create new variable based on survival_check
NP$Time_to_death_or_last_followup_years <- NP$Time_to_survival_check

#Time to ventilation will be considered death (status = dead)
NP$Time_to_death_or_last_followup_years[VentAsDeath] <- NP$Time_to_ventilation[VentAsDeath]
NP$Survival_status[VentAsDeath] <- "dead"

status_no_time<- which( is.na(NP$Time_to_death_or_last_followup_years) & !is.na(NP$Survival_status)) 
time_no_status<- which(! is.na(NP$Time_to_death_or_last_followup_years) & is.na(NP$Survival_status))

#Consider those with survival info but no status as censored at this point
#NP[time_no_status,] #2 people with ALS 1 with ALS/FTD
NP$Survival_status[time_no_status] <- "alive"

#Quickly investigate characteristics of those with status but no survival time
#table(NP$Phenotype[status_no_time]) #3000 People with status no time are Control
# has_onset<-which( is.na(NP$Time_to_death_or_last_followup_years) & !is.na(NP$Survival_status) & !is.na(NP$Age_at_onset_days))
# 
# #table(NP$Phenotype[has_onset]) - all are cases
# NP$Age_at_onset_days[has_onset]
# NP$Age_at_survival_check_days[has_onset]
# NP$Age_at_death_days[has_onset]


#------------# 

####Computing: "survival_status_bin"     
NP<-    NP %>%
  mutate(survival_status_bin = case_when(Survival_status == "alive" ~ 0,
                                         Survival_status == "dead" ~ 1
  ))


#Write out the file
write.csv(NP,file=file.path(path,paste0("New_PM_Pheno",Sys.Date(),".csv")),row.names = FALSE)
