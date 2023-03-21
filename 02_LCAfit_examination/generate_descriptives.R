suppressPackageStartupMessages({
  library(writexl) #For writing the summaries to an xlsx workbook
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
  library(survival)
})

option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file for which results will be generated")
)
opt = parse_args(OptionParser(option_list=option_list))

#Obtain directory tree minus model name. This is where output files will be saved
path <- file.path(dirname(opt$resultsfile),"descriptives/")
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import dataout
data <- fread(file=opt$resultsfile,header=TRUE)

######
### Create descriptive table overall and by class
######

#Recode the cluster variable into a factor
data$C <-  as.factor(data$C)

## Prepare the permanent columns first
#Summarise without grouping by C
unfiltered <- data %>% 
  summarise(C        ='Collapsed',
            n        =n(),
            female   =paste0(sum(Sex_at_birth == 2,na.rm = TRUE)," (",round(mean(Sex_at_birth,na.rm=TRUE)-1,2),")"),
            Bulbar   =paste0(sum(Site_of_Onset == 2,na.rm = TRUE)," (",round(mean(Site_of_Onset,na.rm=TRUE)-1,2),")"),
            mean_sd_ageons   =paste0(round(mean(Age_at_onset_years,na.rm=TRUE),2)," (",round(sd(Age_at_onset_years,na.rm=TRUE),2),")")
  )

#Create summary table filtered by C
filtered <- data %>% 
  group_by(C) %>% 
  summarise(n        =n(),
            female   =paste0(sum(Sex_at_birth == 2,na.rm = TRUE)," (",round(mean(Sex_at_birth,na.rm=TRUE)-1,2),")"),
            Bulbar   =paste0(sum(Site_of_Onset == 2,na.rm = TRUE)," (",round(mean(Site_of_Onset,na.rm=TRUE)-1,2),")"),
            mean_sd_ageons   =paste0(round(mean(Age_at_onset_years,na.rm=TRUE),2)," (",round(sd(Age_at_onset_years,na.rm=TRUE),2),")")
  )

#Combine
sumdat <- rbind(unfiltered,filtered)

## Repeat for: diagnostic delay
if(!is.null(data$Diagnostic_delay_years)){
  
  unfiltered <- data %>%  #Summarise without grouping by C
    summarise(
      med_rng_delay    =paste0(round(median(Diagnostic_delay_years,na.rm=TRUE),2)," (",round(min(Diagnostic_delay_years,na.rm=TRUE),2),", ",round(max(Diagnostic_delay_years,na.rm=TRUE),2),")"),
      med_iqr_delay    =paste0(round(median(Diagnostic_delay_years,na.rm=TRUE),2)," (",round(quantile(Diagnostic_delay_years,na.rm=TRUE)[2],2),", ",round(quantile(Diagnostic_delay_years,na.rm=TRUE)[4],2),")")
    )
  
  filtered <- data %>%  #Summarise grouped by C
    group_by(C) %>% 
    summarise(
      med_rng_delay    =paste0(round(median(Diagnostic_delay_years,na.rm=TRUE),2)," (",round(min(Diagnostic_delay_years,na.rm=TRUE),2),", ",round(max(Diagnostic_delay_years,na.rm=TRUE),2),")"),
      med_iqr_delay    =paste0(round(median(Diagnostic_delay_years,na.rm=TRUE),2)," (",round(quantile(Diagnostic_delay_years,na.rm=TRUE)[2],2),", ",round(quantile(Diagnostic_delay_years,na.rm=TRUE)[4],2),")")
    ) %>%
    dplyr::select(-C)
  
  #Append
  sumdat <- rbind(unfiltered,filtered) %>%
    cbind(sumdat,.)
  
}

## Repeat for: per-country standardised diagnostic delay
if(!is.null(data$CTYDEL)){

  unfiltered <- data %>%  #Summarise without grouping by C
    summarise(
      med_rng_ctydel    =paste0(round(median(CTYDEL,na.rm=TRUE),2)," (",round(min(CTYDEL,na.rm=TRUE),2),", ",round(max(CTYDEL,na.rm=TRUE),2),")"),
      med_iqr_ctydel    =paste0(round(median(CTYDEL,na.rm=TRUE),2)," (",round(quantile(CTYDEL,na.rm=TRUE)[2],2),", ",round(quantile(CTYDEL,na.rm=TRUE)[4],2),")")
    )
  
  filtered <- data %>%  #Summarise grouped by C
    group_by(C) %>% 
    summarise(
      med_rng_ctydel    =paste0(round(median(CTYDEL,na.rm=TRUE),2)," (",round(min(CTYDEL,na.rm=TRUE),2),", ",round(max(CTYDEL,na.rm=TRUE),2),")"),
      med_iqr_ctydel    =paste0(round(median(CTYDEL,na.rm=TRUE),2)," (",round(quantile(CTYDEL,na.rm=TRUE)[2],2),", ",round(quantile(CTYDEL,na.rm=TRUE)[4],2),")")
    ) %>%
    dplyr::select(-C)
  
  #Append
  sumdat <- rbind(unfiltered,filtered) %>%
    cbind(sumdat,.)
}

## Repeat for: disease duration
# Additionally create a survfit estimate of disease duration
if(!is.null(data$Time_to_death_or_last_followup_years)){
  
  unfiltered <- data %>%  #Summarise without grouping by C
    summarise(
      ncensored        =paste0(sum(survival_status_bin == 0,na.rm = TRUE)," (",round(1-mean(survival_status_bin,na.rm=TRUE),2),")"),
      med_rng_survival =paste0(round(median(Time_to_death_or_last_followup_years,na.rm=TRUE),2)," (",round(min(Time_to_death_or_last_followup_years,na.rm=TRUE),2),", ",round(max(Time_to_death_or_last_followup_years,na.rm=TRUE),2),")"),
      med_iqr_survival =paste0(round(median(Time_to_death_or_last_followup_years,na.rm=TRUE),2)," (",round(quantile(Time_to_death_or_last_followup_years,na.rm=TRUE)[2],2),", ",round(quantile(Time_to_death_or_last_followup_years,na.rm=TRUE)[4],2),")")
    )
  
  filtered <- data %>%  #Summarise grouped by C
    group_by(C) %>% 
    summarise(
      ncensored        =paste0(sum(survival_status_bin == 0,na.rm = TRUE)," (",round(1-mean(survival_status_bin,na.rm=TRUE),2),")"),
      med_rng_survival =paste0(round(median(Time_to_death_or_last_followup_years,na.rm=TRUE),2)," (",round(min(Time_to_death_or_last_followup_years,na.rm=TRUE),2),", ",round(max(Time_to_death_or_last_followup_years,na.rm=TRUE),2),")"),
      med_iqr_survival =paste0(round(median(Time_to_death_or_last_followup_years,na.rm=TRUE),2)," (",round(quantile(Time_to_death_or_last_followup_years,na.rm=TRUE)[2],2),", ",round(quantile(Time_to_death_or_last_followup_years,na.rm=TRUE)[4],2),")")
    ) %>%
    dplyr::select(-C)
  
  #Append
  sumdat <- rbind(unfiltered,filtered) %>%
    cbind(sumdat,.)
  
  
  ###### Generate a survfit and extract median
  #Per class
  sv<- survfit(Surv(Time_to_death_or_last_followup_years,survival_status_bin) ~ C,data=data)
  sv_survival<- summary(sv)$table[,"median"]
  names(sv_survival) <- gsub("C=","",names(sv_survival))
  
  #For data overall
  sv_collapse<- survfit(Surv(Time_to_death_or_last_followup_years,survival_status_bin) ~ 1 ,data=data)
  sv_coll_survival<- summary(sv_collapse)$table["median"]
  names(sv_coll_survival) <- "Collapsed"
  
  #Combine
  sv_survival<- c(sv_coll_survival,sv_survival)

  #Append survival as estimated within Surv-fit
  sv_survival <- sv_survival %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    rename(C=rowname,Survfit_median=".")
    
  sumdat<- full_join(sumdat,sv_survival,by="C")
}

## Repeat for: Phenotype
if(!is.null(data$Phenotype)){
  unfiltered <- data %>%  #Summarise without grouping by C
    summarise(
      Phenotype1 = paste0(sum(Phenotype == 1,na.rm = TRUE)," (",round(sum(Phenotype == 1,na.rm = TRUE)/n(),2),")"),
      Phenotype2 = paste0(sum(Phenotype == 2,na.rm = TRUE)," (",round(sum(Phenotype == 2,na.rm = TRUE)/n(),2),")"),
      Phenotype3 = paste0(sum(Phenotype == 3,na.rm = TRUE)," (",round(sum(Phenotype == 3,na.rm = TRUE)/n(),2),")"),
          )
  
  filtered <- data %>%  #Summarise grouped by C
    group_by(C) %>% 
    summarise(
      Phenotype1 = paste0(sum(Phenotype == 1,na.rm = TRUE)," (",round(sum(Phenotype == 1,na.rm = TRUE)/n(),2),")"),
      Phenotype2 = paste0(sum(Phenotype == 2,na.rm = TRUE)," (",round(sum(Phenotype == 2,na.rm = TRUE)/n(),2),")"),
      Phenotype3 = paste0(sum(Phenotype == 3,na.rm = TRUE)," (",round(sum(Phenotype == 3,na.rm = TRUE)/n(),2),")"),
    ) %>%
    dplyr::select(-C)
  
  #Append to main summary
  sumdat <- rbind(unfiltered,filtered) %>%
    cbind(sumdat,.)
}

#Write the summaries to a table
writexl::write_xlsx(x=sumdat, path=paste0(path,paste0("descriptives_",Sys.Date(),".xlsx")),
                    col_names = TRUE)


#=================

### Optionally also generate a violin plot for disease duration and survival
if(!is.null(data$Diagnostic_delay_years)){
  if(!is.null(data$Time_to_death_or_last_followup_years)){
    meltdata<-reshape2::melt(data, id.vars=c("new_num_id","C"), measure.vars=c("Time_to_death_or_last_followup_years","Diagnostic_delay_years"))
  } else {
    meltdata <- reshape2::melt(data, id.vars=c("new_num_id","C"), measure.vars=c("Diagnostic_delay_years"))
  }
  
  violin <- ggplot(meltdata, aes(x=as.factor(C),fill=variable,y=value)) + 
    geom_violin(trim=TRUE)+
    theme_bw()+
    xlab("Class")+
    ylab("Time (years)")+
    theme(legend.position="bottom")
  
  #Save the figure
  ggsave(plot=violin, 
         filename = paste0(path,"Violinplot_",Sys.Date(),".pdf"),
         units="mm",width=200,height=150)

}
