#### Apply Fisher's Exact Test to analyse the distribution of rare variants across classes

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(data.table)
  library(epitools)
  library(optparse)
})

option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file for which results will be generated"),
  make_option("--Predictorcols", action="store", default=NA, type='character',
              help="Specify comma separated string of column names indicating the predictor variables")
)

opt = parse_args(OptionParser(option_list=option_list))

#Import the genetic variables
indivtests <- strsplit(opt$Predictorcols,",")[[1]] 

#Obtain directory tree minus model name. This is where output files will be saved
path <- file.path(dirname(opt$resultsfile),"FishersExact")
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import dataset
data <- data.table::fread(file=opt$resultsfile,
                          select=c("C",indivtests),
                          header=TRUE)

data$C <- as.factor(data$C)


## Specify function for main fishers exact test, allowing for recursive implementation to derive empirical P-values
doFishers <- function(data){
  
  #Call the dataframe 'genes'
  genes <- tibble(data) 
  
  #Create contingency tables across each column of genes
  tabAndName<- function(x,y){
    tab <- table(x,y)
    #1 row can occur if there are no variants whatsoever
    if(nrow(tab)==1){
      tab <- rbind(rep(0,length(tab)),t(tab))
    } else {
      #Convert to matrix and invert the row order
      tab <- matrix(tab,nrow=nrow(tab)) 
      tab <- tab[nrow(tab):1,] #invert the row order
    }
    colnames(tab) <- paste0("Class",levels(y))
    rownames(tab) <- c("Variant","NoVariant")
    
    return(tab)
  }
  
  #Extract all genetic columns and crosstabulate with Class
  genelist <- lapply(genes[colnames(genes)!="C"], tabAndName,y=genes[["C"]]) 
  
  ## Loop across contingency tables and obtain 2x2 class vs other matrices
  for(i in 1:length(genelist)){
    contin <- genelist[[i]] 
    sublist <- vector(mode="list",length=ncol(contin))
    for(j in 1:ncol(contin)){
      twobytwo<- matrix(NA_real_,nrow=2,ncol=2)
      twobytwo[,1] <- contin[,j]
      twobytwo[,2] <- rowSums(contin[,-j])
      colnames(twobytwo) <- c(colnames(contin)[j],"Other")
      rownames(twobytwo) <- c("Variant","NoVariant")
      sublist[[j]] <- twobytwo
      names(sublist)[j] <- paste0(names(genelist)[i],"_2by2_", colnames(contin)[j])
    }
    #Append the 2x2 list to the main genelist
    genelist <- c(genelist, sublist)
  }
  #Run Fishers tests across all contingency tables and save the p-values for all tables in a summary column
  #Epitools fishers exact tests do not compute an overall p-value for rx2 grids.
  #Therefore use Rstats implementation for p-values and epitools for 2x2 odds ratios
  Summary <- data.frame(name=names(genelist), pval=rep(NA_real_,length(genelist)))
  for(f in 1:length(genelist)){
    if(ncol(genelist[[f]])>2){
      fishy<- fisher.test(genelist[[f]],hybrid = TRUE)
    } else {
      fishy<- fisher.test(genelist[[f]])
    }
    Summary[f,2] <- fishy$p.value
    
  }
  genelist$Summary <- Summary
  return(genelist)
}

#Apply function to the real dataset
realresult<- doFishers(data)

######
### ATTAIN ODDS RATIOS for 2x2 tables
#######
#Create elements in the summary list element to store odds ratios and CI
realresult$Summary$upper <- realresult$Summary$lower <- realresult$Summary$OR <- NA_real_
realresult$Summary$adjustedOR <- NA #Include Logical column to specify whether an adjustment was made to the contingency table on accounts of empty cells

#pattern match 2x2 comparisons and determine OR
twobytwo_list <- grep("2by2",names(realresult)) 
for(i in twobytwo_list){
  contingency <- realresult[[i]]
  
  #reference conditions are no variant and 'other' class
  orf<- oddsratio.fisher(x=contingency,rev="both")
    realresult$Summary$lower[i] <- orf$measure["Variant","lower"] #Save the CI irrespective of contingency characteristics
    realresult$Summary$upper[i] <- orf$measure["Variant","upper"]
  
  #If there are any empty cells for the Variant group, perform adjustment and recalculate the OR point estimate
  if(any(contingency==0)){ 
    empty<- which(contingency==0,arr.ind = TRUE)
    
    if(empty[1]==2){
      warning("The NoVariant group contains an empty cell for ",names(realresult)[i])
    }
    contingency[2,empty[2]] <- contingency[2,-empty[2]] #Assign denominator from non-empty class to the empty cells
    contingency <- contingency+1 #Increase all cell counts by 1 (integer used )
    
    orf_adj <- oddsratio.fisher(x=as.matrix(contingency),rev="both")
    
    realresult$Summary$OR[i] <- orf_adj$measure["Variant","estimate"] #If there is no need to adjust, take the original OR  
    realresult$Summary$adjustedOR[i] <- TRUE #Flag that an adjustment was made
  } else {
    realresult$Summary$OR[i] <- orf$measure["Variant","estimate"] #If there is no need to adjust, take the original OR
  }
}

######
### Attain various FDR corrected p-values
######

#FDR correct across all comparisons
realresult$Summary$fdr_pval <- p.adjust(realresult$Summary$pval,method="fdr")

#FDR correct across all nclasses x 2 comparisons
realresult$Summary$fdr_global_pval <- NA
realresult$Summary$fdr_global_pval[-twobytwo_list] <- p.adjust(realresult$Summary$pval[-twobytwo_list],method="fdr")

#FDR correct across all 2 x 2 comparisons
realresult$Summary$fdr_2by2global_pval <- NA
realresult$Summary$fdr_2by2global_pval[twobytwo_list] <- p.adjust(realresult$Summary$pval[twobytwo_list],method="fdr")

#FDR correct family-wise across all 2 x 2 comparisons
#c.f. https://stackoverflow.com/questions/53101814/batch-wise-for-loops-in-r
realresult$Summary$fdr_2by2family_pval <- NA
  adjustbatch<-nrow(realresult$Summary[-twobytwo_list,]) #identify number of global comparisons made. Batches are this -1
  batchsize<- length(levels(data$C)) #Identify the size of each batch according to the number of classes in each family
for (i in 0:(adjustbatch-1)){
  rowBatch <- (batchsize*i + 1):(batchsize * (i + 1))+adjustbatch ##  define the rownumbers of the current batch
  
  realresult$Summary$fdr_2by2family_pval[rowBatch] <- p.adjust(realresult$Summary$pval[rowBatch],method="fdr")
}

######
### Generate qqplot for p-values, based upon: #https://gist.github.com/slowkow/904157
######
gg_qqplot <- function(ps, ci = 0.95) {
  ncomparisons<- length(ps)
  
  qqreal <- data.frame(
    observed=-log10(sort(ps)),
    expected=-log10(ppoints(ncomparisons)),
    clower   = -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:ncomparisons, shape2 = ncomparisons:1)),
    cupper   = -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:ncomparisons, shape2 = ncomparisons:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggqq <- ggplot(qqreal,aes(x=expected,y=observed))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, alpha = 0.5)+
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    theme_bw()+
    xlab(log10Pe) +
    ylab(log10Po)
  
  return(ggqq)
}
#Can be compared to:
#gap::qqunif(realresult$Summary$pval)

qqrealplot<- gg_qqplot(realresult$Summary$pval)

#Save the figure
ggsave(plot=qqrealplot, 
       filename = file.path(path,paste0("qq_real_pval_",Sys.Date(),".pdf")),
       units="mm",width=150,height=150)

######
### Using randomly permuted data, estimate Empirical P across 1000 shuffles of C
######

shuffle <- function(x,data){
  pdata <- data #Take the dataset
  pdata$C <- sample(pdata$C) #Reassign classes randomly
  return(pdata) #output shuffled data
}
#Shuffle data as per function, and save in permute list
permute <- vector(mode="list",length=1000)
permute <- lapply(permute,shuffle,data) 

perm_res<- lapply(permute,doFishers) #apply the fisher function

#Take the p-values from each permuted sample, storing in the pmatrix object
pmatrix <- matrix(NA_real_,nrow=length(perm_res), ncol=nrow(realresult$Summary))
colnames(pmatrix) <- realresult$Summary$name
for(i in 1:length(perm_res)){
  pmatrix[i,]  <- perm_res[[i]]$Summary$pval 
}

#Functions to sum the number of significant values and the proportion significant
colsig <- function(x){sum(x<0.05)}
colsig.prop <- function(x){sum(x<0.05)/length(x)}

summary <- data.frame(pmatrix) %>%
  rename_with(~gsub("_",".",.)) %>% #Remove any underscores from names since these are used for splitting downstream
  summarise(
    across(everything(),list(meanval=mean,sd=sd,median=median,
                             q1=~quantile(.,probs=0.25),q3=~quantile(.,probs=0.75), min=min,
                             max=max,nsig05=colsig,propsig05=colsig.prop))
  ) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_sep = "_",
                      names_to  = c("name", ".value")) %>%
  mutate(name=gsub("\\.","_",name)) #Reinsert the underscores


#Compare the real results with the empirical P summary results
nsmallvsreal <- empirical_P <- vector("numeric",length=nrow(realresult$Summary))
for (i in 1:nrow(realresult$Summary)){
  nsmallvsreal[i] <-  sum(pmatrix[,i]<realresult$Summary$pval[i])
  empirical_P[i] <-  sum(pmatrix[,i]<realresult$Summary$pval[i])/nrow(pmatrix)
}
summary$nsmallvsreal <- nsmallvsreal #Number where the Empirical was smaller than the real result
summary$empirical_P <- empirical_P   #Proportion where the Empirical was smaller than the real result

#Write summary of shuffling to file
data.table::fwrite(x=summary,
                   file=file.path(path,paste0("Simulation_summary_",Sys.Date(),".csv")),
                   col.names = TRUE)

#Combine Empirical P-values with the main comparison list
realresult$Summary <- full_join(realresult$Summary,summary[c("name","empirical_P")],by="name")

#Write fisher's test results and the full summary table
writexl::write_xlsx(x=lapply(realresult,as.data.frame),
                    path=file.path(path,paste0("fishers_",Sys.Date(),".xlsx")),
                    col_names = TRUE)


######
### Lastly, produce tidy summary of main p-values and results 
######

grouptests <- paste0(indivtests,"_2by2")
grouptests_index <- lapply(grouptests,grep,realresult$Summary$name)
nclasses<- length(grouptests_index[[1]])

store<- rowsummary <- vector(mode="numeric")
for(i in 1:length(indivtests)){
  
  indiv_index<- which(realresult$Summary$name==indivtests[i]) #index the matching pvalue
  
  rowsummary[1:5] <- realresult$Summary[indiv_index,c("name","pval","fdr_pval","fdr_global_pval","empirical_P")]

  #index group list for the individual element and extract family-wise p-values
  rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- realresult$Summary$fdr_2by2family_pval[grouptests_index[[indiv_index]]]
  
  #index group list for the individual element and extract fdr-corrected vals across all 2x2 comparisons
  rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- realresult$Summary$fdr_2by2global_pval[grouptests_index[[indiv_index]]]
  
  #index group list for the individual element and extract fdr-corrected vals across ANY comparison
  rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- realresult$Summary$fdr_pval[grouptests_index[[indiv_index]]]
  
  rowsummary[(length(rowsummary)+1):(length(rowsummary)+nclasses)] <- paste0(round(realresult$Summary$OR[grouptests_index[[indiv_index]]],2)," [",
                                                                         round(realresult$Summary$lower[grouptests_index[[indiv_index]]],2), ", ",
                                                                         round(realresult$Summary$upper[grouptests_index[[indiv_index]]],2), "]")
  
  store<- rbind.data.frame(store, rowsummary,stringsAsFactors = FALSE)
  
  rowsummary <- vector(mode="numeric") #reset the length of rowsummary so that the loop doesnt continue appending elements
  
}

colnames(store) <- c("GeneticEff","pval","fdr_pval","fdr_nby2_pval","EmpiricalP",
                     paste0("family_fdr_Class",1:nclasses),
                     paste0("2by2_fdr_Class",1:nclasses),
                     paste0("any_fdr_P_Class",1:nclasses),
                     paste0("OR_Class",1:nclasses))


#Write the summaries to a table
data.table::fwrite(x=store,
                   file=file.path(path,paste0("Final_summary_tab_",Sys.Date(),".tsv")),
                   sep="\t",
                   col.names = TRUE)



##########

