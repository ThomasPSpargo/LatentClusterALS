######
## Combine unique STRENGTH and Project MinE samples to determine a per-country standardised diagnostic delay reference panel
#####
library(dplyr)

#Specify path for input and output files
path <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets"


## Custom Function to pool across country-wise mean and sd and attain a 'total' standardised score
#df is a dataframe with the columns N, mean, and sd
averageMeanSD<- function(df){
  #Extract columns
  mean1 <- df$mean
  sd1 <- df$sd
  n1<- df$N
  
  #Calculate weighted mean for each group and aggregate SD
  wMean1 <- mean1*n1
  wSd1 <-sd1^2*(n1-1)+((wMean1^2)/n1)
  
  #get the total n, sum of weighted means, grouped x2
  totn <- sum(n1)
  sumwMean1 <- sum(wMean1)
  sumwSd1 <- sum(wSd1)
  
  end<- list(N=totn,
             mean=sumwMean1/totn,
             sd= sqrt((sumwSd1-sumwMean1^2/totn)/(totn-1))
  )
  return(end)
}
######
### Derive total average for each dataset and then pool
######

PM_panel<- read.table(file.path(path,"PM_diagDelay_country_standardised.tsv"),
                      header = TRUE, sep="\t")

STR_panel<- read.table(file.path(path,"STR_diagDelay_country_standardised.tsv"),
                       header = TRUE, sep="\t")

#Generate a pooled 'total' diagnostic delay for each panel
STR_pooled<- tibble(rbind(
  cbind(iso2c=NA,data.frame(averageMeanSD(STR_panel)),ID="Total"),
  STR_panel))

PM_pooled<- tibble(rbind(
  cbind(iso2c=NA,data.frame(averageMeanSD(PM_panel)),ID="Total"),
  PM_panel))


#Harmonise across the two panels 
jointPanel<- PM_pooled
for(i in 1:nrow(STR_pooled)){
  match<- which(jointPanel$ID==STR_pooled$ID[i])
  
  df <- rbind(STR_pooled[i,], jointPanel[match,])
  
  combined <- averageMeanSD(df)
  
  jointPanel[match,"mean"] <- combined["mean"]
  jointPanel[match,"sd"] <- combined["sd"]
  jointPanel[match,"N"] <- combined["N"]
  
}

write.table(jointPanel,file=file.path(path,"diagDelay_panel.tsv"),row.names=FALSE, col.names = TRUE,na=".",sep="\t",quote = FALSE)

