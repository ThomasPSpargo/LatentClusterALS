##### Set up ROC functions, 3 functions are defined here
#The first two functions are called within gridROC to generate a grid of ROC curves
#based on observations and predictions in a multiclass classification problem returned by Caret

#pairwiseROCperformance

#oneVsRestROCperformance

#gridROC




#Dependencies are 
#tidyverse
#pROC




#Function to compute pairwise aucs based on a multiclass.roc calculation
pairwiseROCperformance <- function(preds){
  #preds<- finalPred #For testing
  
  #Generates Multi-class auc, with pairwise-aucs. Hand and Till (2001) definition: https://doi.org/10.1023/A:1010920819831
  multi<- multiclass.roc(response=preds$obs,
                         predictor=as.matrix(preds[grepl(paste0("^",paste0(levels(preds$obs),collapse="|"),"$"), colnames(preds))]))
  multiAUC <- as.numeric(auc(multi)) #Total average AUC
  
  #Create extremely long sensitivities and specificities matrix, resembling a melted data frame
  stats<- data.frame(ControlClass=character(0),CaseClass=character(0),Sensitivity=numeric(0),Specificity=numeric(0))
  aucsumDF <- data.frame(ControlClass=character(0),CaseClass=character(0),AUC=numeric(0))
  for(r in 1:length(multi$rocs)){
    nrow1<- length(multi$rocs[[r]][[1]]$specificities)
    nrow2<- length(multi$rocs[[r]][[2]]$specificities)
    stats_tempDir1 <- data.frame(ControlClass=character(nrow1),CaseClass=character(nrow1),Sensitivity=numeric(nrow1),Specificity=numeric(nrow1))
    stats_tempDir2 <- data.frame(ControlClass=numeric(nrow2),CaseClass=numeric(nrow2),Sensitivity=numeric(nrow2),Specificity=numeric(nrow2))
    
    #Now add the class comparisons    
    stats_tempDir1$Sensitivity <- multi$rocs[[r]][[1]]$sensitivities
    stats_tempDir2$Sensitivity <- multi$rocs[[r]][[2]]$sensitivities
    
    stats_tempDir1$Specificity <- multi$rocs[[r]][[1]]$specificities
    stats_tempDir2$Specificity <- multi$rocs[[r]][[2]]$specificities
    
    #Add the levels
    levs <- strsplit(split = "/", names(multi$rocs)[r])[[1]]
    stats_tempDir1$ControlClass <- levs[1]
    stats_tempDir1$CaseClass <- levs[2]
    aucsumDF[nrow(aucsumDF)+1,1:2] <- levs
    aucsumDF$AUC[nrow(aucsumDF)] <- as.numeric(auc(multi$rocs[[r]][[1]]))
    
    
    stats_tempDir2$ControlClass <- levs[2]
    stats_tempDir2$CaseClass <- levs[1]
    aucsumDF[nrow(aucsumDF)+1,1:2] <- rev(levs)
    aucsumDF$AUC[nrow(aucsumDF)] <- as.numeric(auc(multi$rocs[[r]][[2]]))
    
    stats <- rbind(stats,rbind(stats_tempDir1,stats_tempDir2))
  }
  
  return(list(stats=stats,aucsum=aucsumDF,multiAUC=multiAUC))
  
}

oneVsRestROCperformance <- function(preds){
  
  ####
  # Generate 1 vs other roc comparisons
  ####
  
  #Dummy code true classes. Convert to DF and assign ID column to be dropped after completion
  true_bin <- data.frame(ID=1:length(preds$obs), tr=as.numeric(preds$obs),tr_levels=preds$obs) %>% 
    pivot_wider(names_from = tr_levels, values_from = tr, names_prefix= "tr_",values_fill=0,names_sort=TRUE) %>%
    mutate(across(starts_with("tr_"), ~ case_when(.>=1 ~ 1,
                                                  TRUE ~ 0))) %>%
    dplyr::select(.,-(ID))
  
  #Prepare AUC summary table
  AUCsum <- data.frame(ControlClass=character(0),CaseClass=character(0),AUC=numeric(0))
  
  #Loop across dummy coded classifications and try one vs other
  stats<- data.frame(ControlClass=character(0),CaseClass=character(0),Sensitivity=numeric(0),Specificity=numeric(0))
  for(i in 1:ncol(true_bin)) {
    ClassNUM<- gsub("^tr_","",colnames(true_bin)[i])
    
    rocsum<- roc(response=true_bin[[i]],
                 predictor=preds[[which(grepl(ClassNUM, colnames(preds)))]],
                 levels=c(0,1),
                 direction="<") 
    
    stats <- rbind(stats,
                   data.frame(ControlClass=ClassNUM,CaseClass=ClassNUM,
                              Sensitivity=rocsum$sensitivities,
                              Specificity=rocsum$specificities))
    
    AUCsum <- rbind(AUCsum,
                    data.frame(ControlClass=ClassNUM,
                               CaseClass=ClassNUM,
                               AUC=rocsum$auc)) #extract auc from averaged ROC
    
  }
  return(list(stats=stats,aucsum=AUCsum))
}

####
# Generate pairwise roc panel plot, with onevsRest comparisons along the diagonal
####
## finalPred is either a data frame containing the column 'obs' denoting observed groups,
#   and columns which correspond to predicted probabilities of being in each group
#   OR, this is a list containing several dataframes formatted in this way.
# incCV is a logical value. Indicate TRUE to return global performance plus that of subsamples, indicated in the column "resample" for the finalPreds dataframe. Note that this is not used when supplying a list of dataframes
# diagAlpha is a numeric between 0 and 1 to determine the opacity of shading along the one-group vs all other groups ROC curves along the grid diagonal
# labels is a named character vector where names correspond to the group levels recorded in the data and values are labels to use for the facet panels
# plotAUCvals is a logical which if TRUE will return area under the curve for each plotted curve. Set false to disable this behaviour, which may cause issues if plotting a large number of groups
gridROC <- function(finalPred,incCV=FALSE,diagAlpha=0.2,labels=NULL,plotAUCvals=TRUE){
  
  if("list" %in% class(finalPred)){
    #Loop pairwise and onevs rest comparisons across list elements and concatenate
    listPerform_onevrest<- lapply(finalPred, oneVsRestROCperformance)
    listPerform_pairwise<- lapply(finalPred, pairwiseROCperformance)

    ## Summarise all into data.frames
    summaryCols<- c(colnames(listPerform_onevrest[[1]]$stats),"group")
    template = data.frame(matrix(nrow = 0, ncol = length(summaryCols)))
    colnames(template) <- summaryCols

    summaryCols<- c(colnames(listPerform_onevrest[[1]]$aucsum),"group")
    template2 = data.frame(matrix(nrow = 0, ncol = length(summaryCols)))
    colnames(template2) <- summaryCols

    global_pairwiseROC <- global_oneVRest <- list(stats=template
                                                    ,aucsum=template2)
    for(i in 1:length(finalPred)){
      global_oneVRest$stats <- rbind(global_oneVRest$stats,cbind(listPerform_onevrest[[i]]$stats,group=names(listPerform_onevrest)[i]))
      global_oneVRest$aucsum <- rbind(global_oneVRest$aucsum,cbind(listPerform_onevrest[[i]]$aucsum,group=names(listPerform_onevrest)[i]))

      global_pairwiseROC$stats <- rbind(global_pairwiseROC$stats,cbind(listPerform_pairwise[[i]]$stats,group=names(listPerform_pairwise)[i]))
      global_pairwiseROC$aucsum <- rbind(global_pairwiseROC$aucsum,cbind(listPerform_pairwise[[i]]$aucsum,group=names(listPerform_pairwise)[i]))
    }

    if(incCV){
      message("Inclusion of resamples on grid is not supported when supplying a list object, setting incCV=FALSE")
      incCV <- FALSE #Force CV to false for list inputs
    }
  } else {
    
    ## Generate one-vs-rest rocs [which will go along the diagonal]
    global_oneVRest<- oneVsRestROCperformance(finalPred)
    
    ## Obtain stats for the pairwise ROC globally [to plot on the upper and lower tri]
    global_pairwiseROC<- pairwiseROCperformance(finalPred)
    
    if(incCV){ #Obtain stats for each cross-fold resample
      
      #OnevsRest
      cv_oneVRestROCstats<- lapply(unique(finalPred$Resample),function(x,finalPred){ 
        x_cv <-finalPred[finalPred$Resample==x,]
        oneVsRestROCperformance(x_cv)$stats
      },finalPred)
      
      #Pairwise
      cv_pairwiseROCstats<- lapply(unique(finalPred$Resample),function(x,finalPred){ 
        x_cv <-finalPred[finalPred$Resample==x,]
        pairwiseROCperformance(x_cv)$stats
      },finalPred)
      
      #Combine the one vs rest and pairwise stats
      cvStats<- c(cv_oneVRestROCstats,cv_pairwiseROCstats)
    }
  }
  #Concatenate the results and the auc summaries
  globalStats<- list(global_oneVRest$stats,global_pairwiseROC$stats)
  globalAucSum<- as.data.frame(rbind(global_oneVRest$aucsum,global_pairwiseROC$aucsum))
  
  if(!is.null(labels)){
    #Check the labels make sense for the data
    if(length(intersect(global_oneVRest$stats$ControlClass,names(labels))) != length(labels)) {
      stop("names for labels provided do not match the identified in the data.")
    }
    
    ## Convert levels to useful facet names using labels
    globalStats <- lapply(globalStats,function(x,labels){
      x %>%
        mutate(across(ends_with("Class"),~recode(.,!!!labels)))
    },labels)
    
    if(incCV){ #Repeat for resamples
      cvStats <- lapply(cvStats,function(x,labels){
        x %>%
          mutate(across(ends_with("Class"),~recode(.,!!!labels)))
      },labels)
    }
    
    #And apply to summary
    globalAucSum<- globalAucSum %>%
      mutate(across(ends_with("Class"),~recode(.,!!!labels)))
  }

  #Add a geom_rect to shade the background along the diagonal so as to emphasise that the these are slightly different ROC curves
  panelBG<- expand.grid(ControlClass=unique(globalAucSum$CaseClass),CaseClass=unique(globalAucSum$CaseClass)) %>%
    mutate(background=if_else(ControlClass==CaseClass,1,0))
  
  #Supply string-based aesthetic elements via list
  aesthetics<- list(x="Sensitivity",y="Specificity",
                    colour={if("group" %in% names(globalStats[[1]])) "group" else NULL}) %>%
    lapply(., function(x) if(!is.null(x)){sym(x)})
  
  #Define roc plot, lopping across cv and global stats lists
  gridROC<- ggplot()+
     geom_rect(data=panelBG,aes(alpha=background),xmin = -Inf,xmax = Inf, 
               ymin = -Inf,ymax = Inf,fill="azure4") + ### Add shading layer to diagonal panels using alpha mapping
    scale_alpha(range=c(0,diagAlpha),guide="none")+
    {if(incCV) lapply(cvStats,geom_path,mapping=aes(!!!aesthetics),colour='grey',alpha=0.3)}+ #Add the cross-validation performance curves
    lapply(globalStats,geom_path,mapping=aes(!!!aesthetics),alpha=1)+
    {if(plotAUCvals) 
      ggrepel::geom_text_repel(data=globalAucSum,mapping=aes(x=0.25,y=0.25,label=signif(AUC,3),colour=!!aesthetics$colour),
                               show.legend = FALSE,position="stack",direction="y",box.padding = 0.25)}+ #Add the total performance
    facet_grid(rows=vars(ControlClass),cols=vars(CaseClass))+
    theme_bw()+
    theme(axis.text.x = element_text(angle=315, vjust=.5, hjust=0))+
    scale_x_reverse()+
    geom_abline(intercept=1, slope=1,linetype = "dashed")
  
  #Summarise the total pairwise auc and all the one vs rest aucs
  AUCsum <- rbind(c("Total",global_pairwiseROC$multiAUC),
                  global_oneVRest$aucsum[,c("CaseClass","AUC")])
  
  return(list(AverageAUC=AUCsum,gridROC=gridROC,allAUCs=globalAucSum))
}

