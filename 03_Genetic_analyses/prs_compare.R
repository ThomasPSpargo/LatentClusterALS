#### Apply binary logistic regression models to analyse the relationship between classes and polygenic risk scores (PRS)

suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(broom)
  library(optparse)
})


option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file for which results will be generated"),
  make_option("--Predictorcols", action="store", default=NA, type='character',
              help="Specify comma separated string of column names indicating the PRS predictor variables"),
  make_option("--IncludeNofPCs", action="store", default=0, type='numeric',
              help="Specify number of principal components to include in the model; for Project MinE data,between 0 and 20 and expects columns named 'PC#'"),
  # make_option("--getEmpiricalP", action="store", default=FALSE, type='logical',                                     #Option depreciated but retained for completeness with the commented out script portion
  #             help="Specify whether to compute empirical P values for the tests conducted. Can increase processing time substantially"),
  make_option("--controlData", action="store", default=NA, type='character',
              help="Specify a secondary dataset against which to compare")
  
)

opt = parse_args(OptionParser(option_list=option_list))

### Prepare directory structure

#Flag whether analyses will be regressed upon principle components
if(opt$IncludeNofPCs>0){
  pc <- paste0("PC",1:opt$IncludeNofPCs) #This is the null model predictor terms if regressors are included
  path <- file.path(dirname(opt$resultsfile),paste0("PRS_regPC1_",opt$IncludeNofPCs))
} else {
  path <- file.path(dirname(opt$resultsfile),"PRS_reg")
  pc <- 1 #This is the null model predictor terms if regressors are not
}

#Flag control reference
if(!is.na(opt$controlData)) path <- paste0(path,"_wControl_ref")

#Generate control data
if(!dir.exists(path))dir.create(path,recursive = TRUE)

#Import dataset
data <- data.table::fread(file=opt$resultsfile,header=TRUE)

#Indicate the PRS models for testing
indivtests <- strsplit(opt$Predictorcols,",")[[1]]

#Select C, PRS columns and any PC columns, dropping the column '1' if this dummy was included
modelcols <- c("C",indivtests,pc)
modelcols <- modelcols[!grepl("^1$",modelcols)]

#Import control data and assign a dummy class value of 999, then join
if(!is.na(opt$controlData)){
  controldata <- data.table::fread(file=opt$controlData,
                                   select=modelcols[!grepl("^C$",modelcols)],
                                   header=TRUE)
  controldata$C <- 999
  
  data<- full_join(data,controldata)
}

#Subset columns
dataframe <- tibble(data)[,modelcols]

#Place modelling steps within a function so that they can be applied recursively for Empirical P
doRegression <- function(dataframe){ 
  
  #Name the dataframe 'genes')
  genes <- dataframe
  
  
  #Dummy code the original C column (deleting 'C' in the process); 'row' column is temporarily used to keep each record unique
  genes <- genes %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = C, values_from = C, names_prefix= "C_",names_sort=TRUE) %>%
    mutate(across(starts_with("C_"), ~ case_when(.>=1 ~ 1,
                                                 TRUE ~ 0))) %>%
    dplyr::select(-row)
  
  #If a control group is in the column names, restrict the comparison to class vs control
  if("C_999" %in% colnames(genes)){
    genes <- genes %>%
      mutate(across(starts_with("C_"), ~ case_when(.>=1 ~ 1,
                                                   C_999==1 ~ 0,
                                                   TRUE ~ NA_real_))) %>%
      dplyr::select(-C_999)
  }
  
  #Identify the number of classes across which to run comparisons
  c_levs<- grep("^C_",colnames(genes),value=TRUE)
  
  #Extract the column names for PRS predictor variables
  genecols <- colnames(genes)[which(!(colnames(genes) %in% c(c_levs, grep("PC",colnames(genes),value=TRUE))))]
  
  #Loop across comparisons for each PRS and each class
  null.fit <- glm.fit <- vector(mode="list")
  #compare.fit <- vector(mode="numeric") #Depreciated but retained for reference
  for(g in 1:length(genecols)){
    
    loopname <- paste0(genecols[g],"_",c_levs)
    for(s in 1:length(c_levs)){
      #Extract filtered dataset upon which to run the model (ensures parity to null model)
      mod_dat <- genes[!is.na(genes[[genecols[g]]]),]

      #Fit the main model
      fmla <- as.formula(paste(c_levs[s],paste(c(genecols[g],pc),collapse=" + "),sep = " ~ "))
      glm.fit[[length(glm.fit)+1]] <- glm(fmla,family=binomial,data=mod_dat)
      names(glm.fit)[length(glm.fit)] <- loopname[s]
      
      #Fit the 'null' model which is either no predictors or the PC predictors only
      fmla <- as.formula(paste(c_levs[s],paste(pc,collapse=" + "),sep = " ~ "))
      null.fit[[length(null.fit)+1]] <- glm(fmla,family=binomial,data=mod_dat)
      names(null.fit)[length(null.fit)] <- paste0(c_levs[s],"_step0")
      
      #Depreciated but retained for reference
      # #Compute glm.fit significance relative to null.fit with anova
      # #Column name adjusted for consistency between multinom and binomial glm (multinom depreciated)
      # anv <- anova(null.fit[[length(null.fit)]],glm.fit[[length(glm.fit)]])
      # colnames(anv)<- gsub("Resid. (d|D)f","Resid. Df",colnames(anv))
      # compare.fit[[length(compare.fit)+1]] <- with(anv, pchisq(`Resid. Dev`[1]-`Resid. Dev`[2],`Resid. Df`[1]-`Resid. Df`[2],lower.tail=FALSE)) 
      # names(compare.fit)[length(compare.fit)] <- paste0(loopname[s],"_stattest")
    }
  }
  
  
  ### Generate tidy summaries of the model stats
  tidy.fit <- vector(mode="list")
  for(i in 1:length(glm.fit)){
    
    #Glance fit stats for null.fit model and name to be clear this is a model excludes the PRS only
    null.glance <- broom::glance(null.fit[[i]])
    colnames(null.glance) <- paste0("midmodel.",colnames(null.glance))
    
    tidy.fit<- rbind(tidy.fit,
                     cbind(model=names(glm.fit)[i],    #Model name
                           broom::tidy(glm.fit[[i]]),  #Predictor specific model fit statistics
                           broom::glance(glm.fit[[i]]),#Wider model fit statistics
                           null.glance)
    )
  }
  
  #Generate 95% confidence intervals 
  tidy.fit<- cbind(tidy.fit,
                   Beta.lower=tidy.fit$estimate-1.96*tidy.fit$std.error,
                   Beta.upper=tidy.fit$estimate+1.96*tidy.fit$std.error
  )
  
  return(list(tidy.fit=tidy.fit#,
              #pvalues.sum=compare.fit #Depreciated
              ))
}

#Apply function to the real dataset
result <- doRegression(dataframe)
tidy.fit <- result$tidy.fit

#Drop principal component terms for summary generation
if(opt$IncludeNofPCs>0) tidy.fit <- tidy.fit %>%filter(.,!grepl("PC[0-9]+",term))

modfinal<- file.path(path,paste0("LogistPRS_tidy_",Sys.Date(),".csv"))

#Save tidy fit with model as the first column
write.csv(tidy.fit,file=modfinal)


#### Attain various FDR corrected model p-values (NOTE: coefficient p-values are those used in paper, therefore this section is not run and may contain errors introduced as scripts were updated; retained for reference)
#
# pvalues.sum <- result$pvalues.sum
# 
# #Generate dataframe of P-value correction information,including global fdr correction
# pvalues.sum <- data.frame(model=gsub("_stattest","",names(pvalues.sum)),
#                           p=pvalues.sum,
#                           fdr_pval=p.adjust(pvalues.sum,method="fdr"),
#                           row.names = NULL)
# 
# ind<- 1:length(indivtests)#Index all total comparison rows, which is 1-2by2 rows
# 
# #FDR correct across all nclasses x 2 comparisons
# pvalues.sum$fdr_nby2_pval <- NA
# # pvalues.sum$fdr_nby2_pval[ind] <- p.adjust(pvalues.sum$p[ind],method="fdr")
# 
# #FDR correct across all 2 x 2 comparisons
# pvalues.sum$fdr_2by2_pval <- NA
# # pvalues.sum$fdr_2by2_pval[-ind] <- p.adjust(pvalues.sum$p[-ind],method="fdr")
# 
# #FDR correct family-wise across all 2 x 2 comparisons
# #https://stackoverflow.com/questions/53101814/batch-wise-for-loops-in-r
# 
# pvalues.sum$fdr_2by2family_pval <- NA
# #Correct in batches of the following size:
# batchsize<- length(unique(dataframe$C[!grepl(999,dataframe$C)])) 
# twoby2<- pvalues.sum#[-ind,] #Extract only the 2x2 models to pass through loop
# for (i in 0:(length(ind)-1)){
#   rowBatch <- (batchsize*i + 1):(batchsize * (i + 1)) ##  define the rownumbers of the current batch
#   twoby2$fdr_2by2family_pval[rowBatch] <- p.adjust(twoby2$p[rowBatch],method="fdr") #adjust for these rows
# }
# #pvalues.sum$fdr_2by2family_pval[-ind] <- twoby2$fdr_2by2family_pval #insert results to table
# pvalues.sum$fdr_2by2family_pval <- twoby2$fdr_2by2family_pval #insert results to table
# 
# #Write fisher's test to file
# 
# write.csv(x=pvalues.sum,
#           file=file.path(path,paste0("LogistPRS_pvals_",Sys.Date(),".csv")))
# 
# #############
# 
# 
# #Function to generate qqplot for p-values, based upon:
# #https://gist.github.com/slowkow/904157
# gg_qqplot <- function(ps, ci = 0.95) {
#   ncomparisons<- length(ps)
#   
#   qqreal <- data.frame(
#     observed=-log10(sort(ps)),
#     expected=-log10(ppoints(ncomparisons)),
#     clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:ncomparisons, shape2 = ncomparisons:1)),
#     cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:ncomparisons, shape2 = ncomparisons:1))
#   )
#   log10Pe <- expression(paste("Expected -log"[10], plain(P)))
#   log10Po <- expression(paste("Observed -log"[10], plain(P)))
#   ggqq <- ggplot(qqreal,aes(x=expected,y=observed))+
#     geom_point()+
#     geom_abline(intercept = 0, slope = 1, alpha = 0.5)+
#     geom_line(aes(expected, cupper), linetype = 2) +
#     geom_line(aes(expected, clower), linetype = 2) +
#     theme_bw()+
#     xlab(log10Pe) +
#     ylab(log10Po)
#   
#   return(ggqq)
# }
# #Can be compared to:
# #gap::qqunif()
# 
# #Generate qq plot for the total model pvalues only
# # pqq <- pvalues.sum %>%
# # filter(grepl( "Ctotal", model))
# 
# qqrealplot <- gg_qqplot(pvalues.sum$p)
# 
# #Save the figure
# ggsave(plot=qqrealplot, 
#        filename = file.path(path,"qq_real_pval.pdf"),
#        units="mm",width=150,height=150)
# 
# 
# 
# 
# ####################### Now try permuted data #################
# 
# if(opt$getEmpiricalP==TRUE){
#   shuffle <- function(x,data){
#     pdata <- data #Take the dataset
#     pdata$C <- sample(pdata$C) #Reassign classes randomly
#     return(pdata) #output shuffled data
#     
#   }
#   
#   permute <- vector(mode="list",length=1000) #Generate empty list along which to place the data frames
#   permute<- lapply(permute,shuffle) #Shuffle data as per function, and save in permute list
#   
#   
#   perm_res<- lapply(permute,doRegression) #Recursively fit the model
#   
#   pmatrix <- matrix(NA_real_,nrow=length(perm_res), ncol=nrow(pvalues.sum))
#   colnames(pmatrix) <- pvalues.sum$model
#   for(i in 1:length(perm_res)){
#     pmatrix[i,]  <- perm_res[[i]]$pvalues.sum #Take the p-values from each permuted sample
#   }
#   
#   # # #pvector <- as.vector(pmatrix)
#   # qqsimplot <-  gg_qqplot(as.vector(pmatrix)) #Repeat the ggplot for all simulated comparisons
#   # 
#   # #Save the figure
#   # ggsave(plot=qqsimplot,
#   #        filename = paste0(path,"qq_simulated_pval_",Sys.Date(),".pdf"),
#   #        units="mm",width=150,height=150)
#   
#   #Remove underscores from names since these are used in the pivot longer function
#   colnames(pmatrix) <- gsub("_",".",colnames(pmatrix))
#   
#   pmatrix <- data.frame(pmatrix)
#   
#   #Function to sum the number of significant values
#   colsig <- function(x){
#     sum(x<0.05)
#   }
#   
#   colsig.prop <- function(x){
#     sum(x<0.05)/length(x)
#   }
#   
#   summary <- pmatrix %>%  
#     summarise(
#       across(everything(),list(meanval=mean,sd=sd,median=median,
#                                q1=~quantile(.,probs=0.25),q3=~quantile(.,probs=0.75), min=min,
#                                max=max,nsig05=colsig,propsig05=colsig.prop))
#     ) %>%
#     tidyr::pivot_longer(cols = everything(),
#                         names_sep = "_",
#                         names_to  = c("variable", ".value"))
#   
#   
#   nsmallvsreal <- empirical_P <- vector("numeric",length=nrow(pvalues.sum))
#   for (i in 1:nrow(pvalues.sum)){
#     nsmallvsreal[i] <-  sum(pmatrix[,i]<pvalues.sum$p[i])               #Number smaller than the real p-value
#     empirical_P[i] <-  nsmallvsreal[i]/nrow(pmatrix)                    #P value based on this value
#   }
#   summary$nsmallvsreal <- nsmallvsreal
#   summary$empirical_P <- empirical_P
#   
#   
#   #Write summary of shuffling to file
#   data.table::fwrite(x=summary,
#                      file=file.path(path,paste0("Simulation_summary_",Sys.Date(),".csv")),
#                      col.names = TRUE)
#   
# }  
# 
# ##### NOW Produce tidy summary suitable for usage in the presentation
# grouptests <- paste0(indivtests,"_C_[0-9]+")
# 
# grouptests_index <- lapply(grouptests,grep,pvalues.sum$model)
# nclasses<- length(grouptests_index[[1]])
# 
# store<- rowsummary <- vector(mode="numeric")
# for(i in 1:length(indivtests)){
#   
#   #indiv_index<- which(pvalues.sum$model==indivtests[i]) #index the matching pvalue
#   
#   #rowsummary[1:4] <- pvalues.sum[indiv_index,c("model","p","fdr_pval","fdr_nby2_pval")]
#   rowsummary[1] <- pvalues.sum[i,"model"]
#   
#   #if(opt$getEmpiricalP==TRUE){
#   #  rowsummary[5] <- summary$empirical_P[indiv_index] #insert Empirical p-value
#   #} else {
#   #  rowsummary[5] <- NA
#   #}
#   
#   #index group list for the individual element and extract family-wise p-values
#   rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- pvalues.sum$fdr_2by2family_pval[grouptests_index[[i]]]
#   
#   #index group list for the individual element and extract fdr-corrected vals across all 2x2 comparisons
#   rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- pvalues.sum$fdr_2by2_pval[grouptests_index[[i]]]
#   
#   #index group list for the individual element and extract fdr-corrected vals across ANY comparison
#   rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- pvalues.sum$fdr_pval[grouptests_index[[i]]]
#   
#   #Extract  2 by 2 odds ratios from model fits
#   query <- tidy.fit %>%                        #Query model fit stats
#     filter(.$term!="(Intercept)") %>% #Remove intercept terms from the search
#     .[which(.$model %in% pvalues.sum$model[grouptests_index[[i]]] ),] #Identify which terms are in the 2x2 comparisons being assessed
#   
#   #Paste the fit statistics
#   rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- paste0(round(query$estimate,2)," [",
#                                                                              round(query$Beta.lower,2), ", ",
#                                                                              round(query$Beta.upper,2), "]")
#   
#   store<- rbind.data.frame(store, rowsummary,stringsAsFactors = FALSE)
#   
#   rowsummary <- vector(mode="numeric") #reset the length of rowsummary so that the loop doesnt continue appending elements
#   
# }
# 
# colnames(store) <- c("GeneticEff",#"pval","fdr_pval","fdr_nby2_pval","EmpiricalP",
#                      paste0("family2by2_fdr_Class",1:nclasses),
#                      paste0("global2by2_fdr_Class",1:nclasses),
#                      paste0("any_fdr_P_Class",1:nclasses),
#                      paste0("Beta_95CI",1:nclasses))
# 
# #need to collapse store into a normal data.frame
# 
# #Write the summaries to a table
# data.table::fwrite(x=store,
#                    file=file.path(path,"Final_summary_tab.csv"),
#                    sep=",",
#                    col.names = TRUE)
##########

#Write to console to be passed to the figure generating script
cat(modfinal)
