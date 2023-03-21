#### This script generates a figure to visualise polygenic risk score (PRS) analyses

library(ggplot2)
library(dplyr)

############################################################

#Expect args containing paths to either 2 or 4 model summary files.
#File(s) 1 (and 3) are simple models
#File(s) 2 (and 4) are with covariates
#Files 1 and 2 have the reference group as the other classes
#Files 3 and 4 have the reference group as healthy controls

args <- commandArgs(trailingOnly = TRUE)

#Read in all datasets
data <- lapply(args,read.delim,sep=',',header=TRUE)

#Save outputs in the directory containing the first model
path <- dirname(args[1])

#Clean the datasets
cleandat <- function(data){
  
  # #Remove row numbering
   data$X <- NULL
   
  #Keep the SbayesR model rows
  data.sbr<- data[grep("SBayesR", data$model),]
  
  #Remove any interecept terms
  data.sbr <- data.sbr[-grep("Intercept", data.sbr$term),]
  
  # Tidy class name
  data.sbr$term<- gsub(".*_C_","", data.sbr$model)
  
  #Convert PRS names to full disease names
  data.sbr$model<- gsub("cc_.*","", data.sbr$model) %>%
    gsub("ALS","Amyotrophic lateral sclerosis", .) %>%
    gsub("AZ","Alzheimer's disease", .) %>%
    gsub("FTD", "Frontotemporal dementia",.) %>%
    gsub("PD","Parkinson's disease", .) %>%
    gsub("SZ","Schizophrenia",.)
  
  return(data.sbr)
  
}

## Clean datasets and append name indicating the model
data.sbr <- lapply(data,cleandat)

data.sbr[[1]]$covars <- "Simple"
data.sbr[[2]]$covars  <- "Adjusting for\nPC 1 to 5"

#Combine the simple and adjusted models for the comparisons across classes
fulldata<- rbind(data.sbr[[1]],data.sbr[[2]]) %>%
  mutate(across(covars,as.factor))

#If including the control sample, prepare and add these too
if(length(data.sbr)>2){
  data.sbr[[3]]$covars <- "Simple"
  data.sbr[[4]]$covars  <- "Adjusting for\nPC 1 to 5"
  
  #Combine the simple and adjusted models for the comparisons to healthy controls
  extra_data <- rbind(data.sbr[[3]],data.sbr[[4]]) %>%
    mutate(across(covars,as.factor))
  
  # NOT RUN; Get the case control analysis sample sizes 
  # lapply(unique(extra_data$model),function(x,extra_data){
  #   table(extra_data$nobs[extra_data$model==x],extra_data$term[extra_data$model==x])},extra_data)
  
  #Annotate the data configurations
  fulldata$design <- "Other classes"
  extra_data$design <- "Healthy controls"
  
  #Combine all data
  fulldata <- rbind(fulldata,extra_data) %>%
    mutate(across(design,as.factor))
  
  #Set reference level
  fulldata$design <- relevel(fulldata$design,ref="Other classes")
  
}

#Set reference level
fulldata$covars <- relevel(fulldata$covars,ref="Simple")

#Annotate p-value thresholds
fulldata <- fulldata %>%
  mutate(lab=case_when(p.value<0.0001 ~ "****",
                         p.value<0.001 ~ "***",
                         p.value<0.01 ~ "**",
                         p.value<0.05 ~ "*",
                         p.value<0.1 ~ "+",
                         TRUE ~ ''))
           
#Generate the ggplot
total_plot<- ggplot(fulldata,aes(x=term,y=exp(estimate),color=design, shape=covars))+
  geom_point(position=position_dodge(width = 0.6))+
  geom_errorbar(aes(ymin=exp(Beta.lower),ymax=exp(Beta.upper)),width=0,
                position=position_dodge(width = 0.6))+
  facet_wrap(~model,ncol=3)+
  geom_text(aes(y=0.05,label=ifelse(covars=="Simple" & design=="Other classes",lab,"")),show.legend = FALSE)+ #Only show significance for the simple models
  geom_text(aes(y=0.20,label=ifelse(covars=="Simple" & design=="Healthy controls",lab,"")),show.legend = FALSE)+
  theme_bw()+
  geom_hline(yintercept = 1,lty=2)+ #NULL association line
  ylab("Odds ratio [95% confidence interval]")+
  xlab("Class")+
  scale_color_manual(name='Reference group', values=c(scales::muted('blue'), 'firebrick'))+
  scale_shape_discrete(name='Model')+
  labs(caption = "P: **** < 0.0001 < *** < 0.001 < ** < 0.01 < * < 0.05 < + < 0.1")+
  theme(legend.position = c(1, 0),
        legend.justification = c(1.60,-0.10),
        plot.caption = element_text(hjust = 1,size=10),
  )



##### SAVE THE PLOT
ggsave(file.path(path,"PRS_summaryFig.pdf"),total_plot,device="pdf",units="mm",width=250,height=150)


### Produce summary table for glanced fit statistics, based on main model and null model
# #Not run for final model and not thoroughly tested; may contain issues
# #
# #Nagelkerke R2 determined based on rms::lrm function formulae
# #And eqns from here: https://stats.stackexchange.com/questions/8511/how-to-calculate-pseudo-r2-from-rs-logistic-regression
# 
#   sumtab <- fulldata %>%
#     mutate(chisq = midmodel.deviance-deviance,
#            mod.pval = 1-pchisq(chisq,midmodel.df.residual-df.residual), ##### DOUBLE CHECK THIS
#            pthresh = case_when(mod.pval<0.0001 ~ "****",
#                                mod.pval<0.001 ~ "***",
#                                mod.pval<0.01 ~ "**",
#                                mod.pval<0.05 ~ "*",
#                                mod.pval<0.1 ~ "+",
#                                TRUE ~ ''),
#            nagel_r2 =  (1-exp(-(midmodel.deviance-deviance)/nobs))/(1-exp(2*midmodel.logLik/nobs)),
#            across(c(mod.pval,nagel_r2), signif,digits=2),
#            across(c(chisq), round,digits=2)
#     ) %>%
#     select(covars,model,term,nobs,chisq,mod.pval,pthresh,nagel_r2) %>%
#     #mutate(across(c(chisq,mod.pval,nagel_r2), signif,digits=2)) %>%
#     rename(c(PRS=model,Class=term))
#   
#   tempbind<- sumtab %>%
#     filter(covars!="Simple") %>%
#     select(-c(covars,nobs)) %>%
#     rename(c(mvar.chisq=chisq,mvar.mod.pval=mod.pval),mvar.pthresh=pthresh,mvar.nagel_r2=nagel_r2)
#   
#   sumtab <- sumtab %>%
#     filter(covars=="Simple") %>%
#     select(-covars) %>%
#     left_join(tempbind,by=c("PRS","Class"))
#   
#   write.table(sumtab,file.path(path,"prs_summary_table.csv"),sep=",",row.names = FALSE)
  
