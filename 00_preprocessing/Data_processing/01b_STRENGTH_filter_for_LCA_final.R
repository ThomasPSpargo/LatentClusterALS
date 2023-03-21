library(dplyr) #Required for piping
library(tidyr) #For pivot_wider
library(lubridate) #Manipulate time


#Specify path containing input files
inpath <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/datasets"

#Specify path out for files
outpath <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets"
if(!dir.exists(outpath)){dir.create(outpath)}

######
### Read in files
######
#STRENGTH data
STRdat<- read.csv(file.path(inpath,"combined_95IDs_20150227pp_AI.csv"),
         header=TRUE,stringsAsFactors = FALSE,na.strings = "-999")

#Read in IDs of people included in project mine sample Project MinE [an output file of script 1a]
PMid <- read.table(file.path(outpath,"PM_identities.tsv"),
                 header=TRUE,sep="\t")

######
### Define custom function for id matching with Project MinE
######

##Function to search for pattern matches across a list of vectors or a single vector, Returns TRUE if any match is found.
#pattern is the pattern to find, xmatch is the list or vector against which to compare
idSearch <- function(pattern,xmatch){
  
  if("list" %in% class(xmatch)) {
    xmatch <- unlist(xmatch)          #Convert to a vector across which to search  
  }
  
  anyPresent<- pattern %in% xmatch #Search for match of pattern in all elements of x and return true or false for each element of x
  drop<- any(anyPresent)          #Return TRUE if a match is found with any element
  return(drop)
}


#####
### Recode continuous Strength data
#####

#colnames(STRdat)
#[1] "id_orig"    "id"         "origin"     "sex"        "diagnosis"  "bulbar"     "ethnicity"  "eec"        "short_eec" 
#[10] "fh"         "date_on"    "date_diag"  "diag_delay" "date_death" "last_seen"  "duration"   "status"     "age_onset" 

#Convert delay (months) into delay (years)
STRdat$DELAY <- dmonths(STRdat$diag_delay) %>%
  time_length(.,unit="year")

#Derive time from onset until death or censoring or death based on delay (years) and duration (years).
  #Duration covers time from diagnosis until death
STRdat$SRV <- STRdat$DELAY+STRdat$duration

#Identify people with delay but no duration information (and vice versa), and assign the available variable as their as survival time and then treat them all as censored since some portion of their disease duration is absent
delayNOduration <- which(is.na(STRdat$SRV) & !is.na(STRdat$DELAY) & is.na(STRdat$duration))
durationNOdelay <- which(is.na(STRdat$SRV) & is.na(STRdat$DELAY) & !is.na(STRdat$duration))

STRdat$SRV[delayNOduration] <- STRdat$DELAY[delayNOduration]
STRdat$SRV[durationNOdelay] <- STRdat$duration[durationNOdelay]
STRdat$status[c(delayNOduration,durationNOdelay)] <- 0


####Finally, check for mismatched missing-ness in SRV and SRVBIN,
#code any unpaired records as NA if SRV is missing or censored if SRVBIN is missing, else they will be removed from LCA
STRdat$status[which(!is.na(STRdat$SRV) & is.na(STRdat$status))] <- 0
STRdat$status[which(is.na(STRdat$SRV) & !is.na(STRdat$status))] <- NA

#For people with 0 days from diagnosis to survival/last_seen, add one day of survival
#(required to be used as survival variable)
STRdat$SRV[which(STRdat$SRV<=0)] <- 1/365.25 #Survival must be at least 1 day

#For diagnostic delay, code impossible values as missing
STRdat$DELAY[which(STRdat$DELAY<=0)] <- NA   #Diagnosis before onset is not possible



#Identify people overlapping with project mine datasets
id_list <- unlist(strsplit(PMid$All_known_IDs,split=","))    #Split the patterns of oldIDs to give a list (which is then 'unlist'-ed). Each element represents an ID present within the Project MinE datasets.
compareID <- sapply(STRdat$id_orig,idSearch,id_list)         #Check which strength IDs are within the Project MinE ID list
#length(which(compareID==TRUE)) #1824

# #Determine country distribution for overlapping samples
# STRdat %>%
#   filter(compareID==TRUE) %>%
#   group_by(origin) %>%
#   summarise(n())

#####
### Filter to people not in Project MinE and WITH some diagnostic label and recode variables for parity to the PM dataset
#####
STRdat <- STRdat %>%
  filter(!is.na(diagnosis) & compareID==FALSE) %>% 
  mutate(origin = case_when(origin == 1 ~ "NL", #NL
                            origin == 2 ~ "GB", #UK
                            origin == 3 ~ "IT", #Italy
                            origin == 4 ~ "BE", #Belgium 
                            origin == 5 ~ "IE"  #Ireland
  ),
  sex = case_when(sex == 0 ~ 1, #Male
                  sex == 1 ~ 2  #Female
  ),
  bulbar = case_when(bulbar == 0 ~ 1, #non-bulbar
                     bulbar == 1 ~ 2, #bulbar
                     
  ),
  diagnosis = case_when(diagnosis == 1 ~ 1, #ALS
                    diagnosis == 2 ~ 3, #PMA
                    diagnosis == 3 ~ 2  #PLS
  
  ),
  new_num_id = as.numeric(paste0(555,row_number()))
  ) #Assign unique numerical id number - important for matching mplus output with the raw file)

#Write file containing all fields from the dataset
write.table(STRdat,file=file.path(outpath,"STR_LCA_Allfields.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)

## Write index file for use in identifying samples
STRdat[c("new_num_id",
             "id_orig",
             "origin")] %>%
  write.table(.,file=file.path(outpath,"STR_identities.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)


######
### Generate a STRENGTH-dataset per-country standardisation panel for diagnostic delay
######
delay_panel<- STRdat %>%
  filter(!is.na(DELAY)) %>% 
  group_by(origin) %>%
  summarise(N=n(),
            mean=mean(DELAY),
            sd=sd(DELAY)
  ) %>%
  mutate(ID=origin)  %>%
  rename(iso2c=origin)

write.table(delay_panel,file=file.path(outpath,"STR_diagDelay_country_standardised.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)