#Script for comparing model fits across mplus output directories
library(MplusAutomation)
library(ggplot2)

#Define path to the initial latent class fits to iterate over
#Loop if multiple paths
runpath<- commandArgs(trailingOnly = TRUE)

#Source helper functions
helperFuncs <- "/scratch/users/k1802739/ALS_LCA/scripts/mplusHelper"
invisible(lapply(list.files(helperFuncs,full.names = TRUE,pattern=".R"),source))

graphIC <- function(runpath){
  bestmodels <- readModels(runpath, recursive=FALSE) #Read in the results of initial model runs
  
  #Identify the class number of classes in each model in ascending order
  nclas <- vector(mode="numeric",length=length(bestmodels))
  for(i in 1:length(bestmodels)){
    nclas[i] <- bestmodels[[i]]$summaries$NLatentClasses
  }
  
  #reorder bestmodels to be ascending according to number of classes
  bestmodels <- bestmodels[order(nclas)]
  
  message("generate fit statistics for models in the ", runpath, "directory:", names(bestmodels))
  
  #Identify the k-1 fit from the previous model #Add k-1 LL value to each model
  for(b in 2:length(bestmodels)){
    bestmodels[[b-1]]$summaries$LL
    bestmodels[[b]]$kmin1LL <- bestmodels[[b-1]]$summaries$LL
  }
  
  #Next tabulate summary statistics
  fitsum<- compareFits(bestmodels)
  
  write.table(fitsum, 
              file=file.path(runpath,paste0("comparefit",Sys.Date(),".csv")),
              sep=",",
              row.names = FALSE,
              col.names = TRUE)
  
  #Return plot comparing the figures
  fitstats <- ggplot(fitsum,aes(x=1:nrow(fitsum)))+
    geom_line(aes(y=AIC,color="AIC"))+
    geom_line(aes(y=BIC,color="BIC"))+
    geom_line(aes(y=aBIC,color="aBIC"))+
    xlab("Model number")+
    ylab("Information criterion value")+
    labs(colour="Statistic")+
    theme_bw()+
    ggtitle("Compare information criterion",subtitle=paste0(names(bestmodels)[1]," : ",names(bestmodels)[length(bestmodels)]))
  
  ggsave(plot=fitstats,file.path(runpath,paste0("compareIC",Sys.Date(),".pdf")),height = 30, width = 20, units="cm")
  
  message("----graph done----")
}

#Apply the function across directories listed in runpath
lapply(runpath,graphIC)

