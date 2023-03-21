suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
  library(pROC)
  library(class)
})

option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file from which a model will be obtained"),
  make_option("--keepcols", action="store", default=NA, type='character',
              help="Indicate columns which will be utilised in the script"),
  make_option("--Predictorcols", action="store", default=NA, type='character',
              help="Indicate columns which will be utilised as model predictors"),
  make_option("--dataset2", action="store", default=NULL, type='character',
              help="Specify a secondary dataset [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

### Prepare for analysis

#Prepare file structure
path <- file.path(dirname(opt$dataset2),"KNN")
if(!dir.exists(path)){dir.create(path, recursive = TRUE)}

#Accuracy function
accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

#Flag columns to use in the analysis
keepcols <- strsplit(opt$keepcols,split=",")[[1]]
predictorcols <- strsplit(opt$Predictorcols,split=",")[[1]]

#Import primary and secondary datasets
traindata <- fread(file=opt$resultsfile,select=keepcols, header=TRUE)
traindata$C <- as.factor(traindata$C)
traindata <- na.omit(traindata)

testdata <- fread(file=opt$dataset2,select=keepcols, header=TRUE)
testdata$C <- factor(testdata$C,levels=levels(traindata$C))
testdata <- na.omit(testdata)

### Specify function to run KNN and then
applyKNN <- function(train_df,test_df,classifier="C",predictorcols,path=NULL){
  
  # ##extract column of train dataset to use as 'cl' argument in knn function.
  train_class <- train_df[[classifier]]
  
  # ##extract equiv column if test dataset to measure the accuracy
  test_class <- test_df[[classifier]]
  
  #Drop classifier as a variable in the training and test data
  train_df<- train_df[,predictorcols,with=FALSE]
  test_df <- test_df[,predictorcols,with=FALSE]
  
  
  #Try KNN for k values between 1 and 20. Try each value 20 times to determine the best in average performance
  performance <- matrix(rep(NA_real_,20*20*2), ncol=2)
  for(i in 1:20){ #Loop across number of neighbours
    for(j in 1:20){ #Repeat across multiple seeds
      row <- ((i-1)*20)+j #The matrix row for the output
      knnest<- knn(train_df,test_df,cl=train_class,k=i) ##run KNN for i classes
      
      compare<- table(test_class,knnest) ##create confusion matrix of real vs predicted class
      performance[row,2] <- accuracy(compare) #examine prediction accuracy
      performance[row,1] <- i #examine prediction accuracy
    }
  }
  colnames(performance) <- c("Neighbours","Accuracy")
  
  #Derive performance and the best performing k-neighbours
  performance<- as_tibble(performance)
  meanperformance<- performance %>%
    group_by(Neighbours) %>%
    dplyr::summarise(mean_accuracy=mean(Accuracy))
  best <- which(meanperformance$mean_accuracy==max(meanperformance$mean_accuracy))[1] #Identify best fit in test data
  
  ##run KNN one further time with set seed for the optimum number of neighbours
  set.seed(50)  
  pr <- knn(train_df,test_df,cl=train_class,k=best) 
  
  ##Generate confusion matrix for best fit
  tab <- table(test_class,pr)
  maketab <- matrix(tab,nrow=nrow(tab))
  colnames(maketab) <- paste0("Pred_",1:nrow(maketab))
  rownames(maketab) <- paste0("True_",1:nrow(maketab))
  
  #Store in output performance per k Neighbours and then the crosstabulation of the best fit
  output<- list(Performance= performance,
                crosstab=data.frame(cbind(rows=rownames(maketab),maketab))
  )
  
  if(!is.null(path)){ #If path specified to function, write out xlsx and ggplot files
    writexl::write_xlsx(x=output, path=file.path(path,paste0("KNN_results_",Sys.Date(),".xlsx")), col_names = TRUE)
    
    #Plot performance by k and save
    plotperf<- ggplot(meanperformance,aes(x=Neighbours,y=mean_accuracy))+
      geom_line()+
      geom_point(data=performance,aes(x=Neighbours,Accuracy),alpha=0.3)+
      geom_vline(xintercept=best,lty=2)+
      theme_bw()+
      scale_x_continuous(breaks = seq(0,max(performance$Neighbours), by = 2))+
      ggplot2::labs(x="Number of neighbours considered",y="Percentage of class assignments correctly predicted")
  
    ggsave(plotperf,filename = file.path(path,"KNN_performance.pdf"),
           units="mm",width=150,height=150)
  }
  
  #Return the predicted classes from the best model and the output list
  return(list(pr=pr,true=test_class,output=output)) 
}

#Apply KNN to external dataset, ext$pr and ext$true will carry forward into AUC calculation
ext <- applyKNN(traindata,testdata,"C",predictorcols,path=path)

######
### Determine area under the curve for each class vs rest
######

#Dummy code predicted and actual classes. Convert to DF and assign ID column to be dropped after completion
pr_bin <- data.frame(ID=1:length(ext$pr), pr=as.numeric(ext$pr)) %>% 
  pivot_wider(names_from = pr, values_from = pr, names_prefix= "pr_",values_fill=0,names_sort=TRUE) %>%
  mutate(across(starts_with("pr_"), ~ case_when(.>=1 ~ 1,
                                                  TRUE ~ 0))) %>%
  dplyr::select(.,-(ID))

true_bin <- data.frame(ID=1:length(ext$true), tr=as.numeric(ext$true)) %>% 
  pivot_wider(names_from = tr, values_from = tr, names_prefix= "tr_",values_fill=0,names_sort=TRUE) %>%
  mutate(across(starts_with("tr_"), ~ case_when(.>=1 ~ 1,
                                                TRUE ~ 0))) %>%
  dplyr::select(.,-(ID))

#Check that all true columns have been represented in the predictions, any absent should be dropped
true_bin <- true_bin[which(gsub("^t","p",colnames(true_bin)) %in% colnames(pr_bin))]

#Produce ROC for each class vs all other classes, save this and the AUC for that fit
rocsum <- list()
AUCsum <- matrix(NA_character_,ncol=2,nrow=ncol(true_bin))
for(i in 1:ncol(true_bin)) {
  roc <- roc(response=true_bin[[i]],
             predictor=pr_bin[[i]],
             levels=c(0,1)) #fit ROC
  rocsum[[i]]<- roc
  
  #Name list element by class and AUC
  colnum<- gsub("tr_","",names(true_bin)[i])
  auc <- round(roc$auc,3)
  names(rocsum)[i] <- paste0(colnum," (AUC = ",auc,")") 
  
  AUCsum[i,] <- c(paste0("Class",colnum), auc) #Save auc in separate matrix
}

#Set palette for ROC figure
cbbPalette <- c("gold", "violetred", "#56B4E9", "#000000", "limegreen", "#D55E00", "#CC79A7", "#F0E442")

roc_plot<- ggroc(rocsum)+
  theme_bw()+
  geom_abline(intercept=1, slope=1,linetype = "dashed")+
  scale_color_manual(name="Class vs other",values = cbbPalette)

#Save the figure
ggsave(plot=roc_plot, filename = file.path(path,paste0("roc_plot_",Sys.Date(),".pdf")),
       units="mm",width=150,height=150)

#Save the AUC values separately
data.table::fwrite(data.table(AUCsum), file=file.path(path,paste0("auc_vals_",Sys.Date(),".txt")),
                     sep="\t",col.names = FALSE)


