
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  #library(randomForest) #Called downstream, as required
  #library(xgboost) #Called downstream
  library(optparse)
  library(pROC)
})

option_list = list(
  make_option("--finalTune", action="store", default=NA, type='character',
              help="Specify the source file from which a model will be obtained"),
  make_option("--tuneFor", action="store", default=NA, type='character',
              help="Specify the parameter to tune for"),
  make_option("--scriptDir", action="store", default=NULL, type='character',
              help="Specify a directory containing scripts to source in order to load functions required to run this script")
)

opt = parse_args(OptionParser(option_list=option_list))

options(echo = FALSE)
cat('Options are:\n')
print(opt)

#Load helper functions for ROC curve generation
source(file.path(opt$scriptDir,"roc_func.R"))

#Prepare the output directory
path <- gsub(".Rds$","",opt$finalTune)
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}
  
# Read in the model
model <- readRDS(opt$finalTune)

#Take forward the predictions
finalPred <- model$pred

#Remove the 'X' prefix from the model
finalPred$obs <- as.factor(gsub("^X","",finalPred$obs))
finalPred$pred <- as.factor(gsub("^X","",finalPred$pred))
colnames(finalPred)[grep("X[0-9]+$",colnames(finalPred))] <- gsub("^X","",grep("X[0-9]+$",colnames(finalPred),value=TRUE))


#####
## Next steps differ by nclasses
#####

## Obtain metrics of per-class sample size and of performance (e.g. precision, specificity and sensitivity), averaged over resamples
if(length(levels(model$trainingData$.outcome))>2){
  
  #Logical for determining whether to process a multiclass or binary model
  isMulticlass <- TRUE
  
  #Obtain sample sizes per group
  classNs<- c("Sample size",table(model$trainingData$.outcome))
  names(classNs) <- c("Metric",gsub("^X","Class: ",names(classNs)[2:length(classNs)]))
  
  confMat<- caret::confusionMatrix(data=finalPred$pred,reference=finalPred$obs)
  performance <- t(confMat$byClass) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    rename(Metric=rowname)
  
  performance<- rbind(classNs,performance)
} else {
  #Logical for determining whether to process a multiclass or binary model
  isMulticlass <- FALSE
  
  #Obtain sample sizes per group
  classNs<- table(model$trainingData$.outcome)
  names(classNs) <- gsub("^X","N_Class: ",names(classNs))
  
  #Set the positive class to the 2nd, given it's the 'non-reference'
  confMat<- caret::confusionMatrix(data=finalPred$pred,reference=finalPred$obs,positive="2")
  
  performance<- t(c(classNs, confMat$byClass))
}
write.table(x=performance,file=file.path(path,paste0("Metrics.tsv")),sep="\t",row.names = FALSE)

###### 
### Plot the roc results 
###### 
if(isMulticlass){
  
  #Set a named vector
  labels= setNames(gsub("X","Class",levels(finalPred$obs)), levels(finalPred$obs)) # create named vector
  
  #Call to the roc helper function directly for plotting the multiclass result on a grid
  aveROC <- gridROC(finalPred,incCV=TRUE,labels=labels)
  
  #Save the faceted ggplot for combination with those from other methods
  saveRDS(aveROC$gridROC,file=file.path(path,paste0("GridRoc.Rds")))
  
  ## Plot the one vs rest total ROC curves lines on a single panel
  
  #Define colour palette
  cbbPalette <- c("gold", "violetred", "#56B4E9", "#000000", "limegreen", "#D55E00", "#CC79A7", "#F0E442")
  
  #Compute onevsRest roc performance using custom function
  global_oneVRest<- oneVsRestROCperformance(finalPred)
  
  #Extract data for plotting
  totPlot<- global_oneVRest$stats
  
  #Assign auc values to each class
  classes <- as.character(global_oneVRest$aucsum[,"CaseClass"])
  names(classes) <- paste0(global_oneVRest$aucsum[,"CaseClass"]," (",signif(global_oneVRest$aucsum[,"AUC"],3),")")
  for(i in 1:length(classes)){
    totPlot$CaseClass[totPlot$CaseClass==classes[i]] <- names(classes)[i]
  }
  
  onePanelRoc<- ggplot(totPlot,aes(x=Specificity,y=Sensitivity,color=CaseClass))+
    geom_path()+ #Add the cross-validation performance curves
    scale_color_manual(name="Class (AUC)",values=cbbPalette)+
    geom_abline(intercept=1, slope=1,linetype = "dashed")+
    scale_x_reverse()+
    theme_bw()+
    labs(x="Specificity",y="Sensitivity")
  
  #Save the overall performance for one vs rest and total for all one-vs-rest
 write.table(x=aveROC$AverageAUC,file=file.path(path,paste0("roc_auc.tsv")),sep="\t",row.names = FALSE)
  
} else {
  
  #Generate a standard ROC curve for a binary model
  
  #Obtain global roc
  binRoc_global<- roc(response=finalPred$obs,
                      predictor=finalPred$`2`,
                      levels=c(1,2),
                      direction="<")
  
  #Obtain roc for each cross-fold resample
  binRoc_CVstats <- lapply(unique(finalPred$Resample),function(x,finalPred){ 
    x_cv <-finalPred[finalPred$Resample==x,]
    
    roc(response=x_cv$obs,
        predictor=x_cv$`2`,
        levels=c(1,2),
        direction="<")
  },finalPred)
  
  
  ## Define ROC plot

  # lopping across cv and global stats lists for 'pairwise' and using only the global stat for 'onepanel'
  # Therefore define a 'root' plot for the recycled elements
  rocRoot<- ggplot()+
    annotate("text",x=0.25,y=0.25,label=signif(binRoc_global$auc,3))+
    theme_bw()+
    labs(x="Specificity",y="Sensitivity")+
    scale_x_reverse()+
    geom_abline(intercept=1, slope=1,linetype = "dashed")
  
  #'Pairwise' auc includes the repeat CVs (use list to mirror output of gridROC)
  aveROC <- list()
  aveROC$gridROC<- rocRoot+
    lapply(binRoc_CVstats,function(x){
      geom_path(as.data.frame(x[c("specificities","sensitivities")]),mapping=aes(x=specificities,y=sensitivities),colour='grey',alpha=0.3)
    })+ #Add the cross-validation performance curves
    geom_path(data=as.data.frame(binRoc_global[c("specificities","sensitivities")]),mapping=aes(x=specificities,y=sensitivities),colour='black',alpha=1) #Add the total performance
  
  #'onepanel' gives just the main result
  onePanelRoc<- rocRoot+
    geom_path(data=as.data.frame(binRoc_global[c("specificities","sensitivities")]),mapping=aes(x=specificities,y=sensitivities),colour='black',alpha=1) #Add the total performance
  
}

# #Save the pairwise and one-vs-rest gridwise ROCs plot
ggsave(plot=aveROC$gridROC,
       filename = file.path(path,paste0("Pairwise_roc_",Sys.Date(),".pdf")),
       units="mm",width=150,height=150)

#Save the one-panel ggplot for combination with those for other data configurations
saveRDS(onePanelRoc,file=file.path(path,paste0("OnePanel_roc.Rds")))


######
# Obtain per-model feature importance metrics
######
nclass <- length(model$finalModel$obsLevels)


if(class(model$finalModel)=="xgb.Booster"){
  #For an xgboost model:
  library(xgboost)
  
  ##Extract overall performance then (for multiclass scenarios) loop across every nclass trees to extract Multiclass per-class feature importance, following the xgb.importance function example documentation
  xgbImp <- cbind(class="Overall",xgb.importance(model=model$finalModel))
  
  if(isMulticlass){
    for(i in 0:(nclass-1)){
      xgbImp <- rbind(xgbImp,
                      cbind(class=paste0("X",i+1),
                            xgb.importance(model=model$finalModel,
                                           trees = seq(from=i,by=nclass,length.out=model$finalModel$niter))
                      )
      )
    }
  }
  #Define the XGboost performance metrics
  metrics <-c("Gain","Cover","Frequency")
  
  #Prepare for plotting
  impSummary<- xgbImp %>%
    pivot_longer(cols=all_of(metrics),names_to = "Metric",values_to = "value") %>%
    mutate(value=value*100) #Rescale as a percentage
  
} else {
  #For a random forest model:
  library(randomForest)
  
  #Define the rf performance metrics [This is necessary for the factor transformation; must be comparable to xgboost options]
  metrics <- c("Random Forest")
  
  # #This calculation gives the same as varImp from caret
  # imp<- randomForest::importance(model$finalModel)[,1:nclass]
  # imp <- imp/max(imp)
  
  #Prepare for plotting; caret extracts and scales the RF model to 100 automatically
  impSummary<- varImp(model)$importance %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    rename(Feature=rowname) %>%
    mutate(Overall=rowMeans(across(starts_with("X")))) %>%
    pivot_longer(cols=!Feature,names_to = "class",values_to = "value") %>%
    mutate(Metric=metrics)
  
  #For binary classification outcomes, filter to only the 'overall' importance since that for X1 and X2 is universal
  #Logic check for identical is implemented to ensure this introduces no unforeseen miscalculations
  if(!isMulticlass && identical(impSummary[impSummary$class=="Overall","value"],impSummary[impSummary$class=="X1","value"])){
    impSummary <- impSummary %>%
      filter(class == "Overall")
  }
  
}

#Declare plot-friendly names for features
recodes=c("Time_to_death_or_last_followup_years"="Disease duration",
"Age_at_onset_years"="Age of onset",
"CTYDEL"="Diagnostic delay",
"Sex_at_birth" = "Sex",
"Site_of_Onset" = "Site of onset",
"Phenotype_1" = "Clinical diagnosis (ALS)",
"Phenotype_2" = "Clinical diagnosis (PLS)",
"Phenotype_3" = "Clinical diagnosis (PMA)",
"c9_status" = "C9orf72 expansion",
"ALScc_SBayesR" = "ALS PRS",
"FTDcc_SBayesR" = "FTD PRS",
"SZcc_SBayesR" = "Schizophrenia PRS",
"PDcc_SBayesR" = "Parkinson's disease PRS",
"AZcc_SBayesR" = "Alzheimer's disease PRS",
"AutoProt" = "Autophagy & proteostasis",
"CytoTransp" = "Cyt. Dynam. & Ax. Transp.",
"RNAFunc" = "RNA function")




#Save summary to file which could be used for generating a feature importance plot across models
write.table(x=impSummary,file=file.path(path,paste0("ImpSummary.tsv")),sep="\t",row.names = FALSE)  

#### Plot variable importance the single model using readable names that 
#Prepare data and pass to ggplot via pipe
Opt_varimpplot <- impSummary %>%
  mutate(Feature=recode(Feature,!!!recodes),
         class = relevel(as.factor(gsub("^X","Class ",class)),ref="Overall"),
         Feature = factor(Feature,levels=Feature[class=="Overall" & Metric==metrics[1]][order(value[class=="Overall" & Metric==metrics[1]])]),
         valuelabel=round(value, 2)
  ) %>%
  ggplot(aes(x = value, y = Feature,fill=Metric))+
  facet_wrap(vars(class))+
  geom_col(position=position_dodge())+
  #geom_text(aes(label = ifelse(value>max(value)*0.40,valuelabel,"")), hjust=1.1, color="white", size=4) +
  #geom_text(aes(label = ifelse(value<=max(value)*0.40,valuelabel,"")), hjust=-0.1, color="black", size=4) +
  theme_bw()+
  labs(x = "Relative Importance",y = "Feature")

if(length(metrics)==1){
  Opt_varimpplot <- Opt_varimpplot+
    scale_fill_manual(values="grey40",guide='none')
}

ggsave(plot=Opt_varimpplot,
       filename = file.path(path,paste0("varImportance_opt",Sys.Date(),".pdf")),
       units="mm",width=150,height=150)



