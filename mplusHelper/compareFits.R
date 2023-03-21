## Helper function to calculate several comparison of fit statistics for a selection of related Mplus models

#http://www.statmodel.com/download/LCAFAQ.pdf
compareFits <- function(models){
  
  fits <- data.frame(model=NA,
                     AIC=NA,
                     BIC=NA,
                     aBIC=NA,
                     sic=rep(NA,length(models))
                     )

  for(i in 1:length(models)){
    mod<- models[[i]]
    
    #Separate this to keep numeric class in cols 2 onwards
    fits[i,1]            <- mod$input$title
    fits[i,2:ncol(fits)] <- c(mod$summaries$AIC,mod$summaries$BIC,mod$summaries$aBIC,-.05*(mod$summaries$BIC))
  }
  

  
  sic <- fits$sic #Take Sic for every model
    BF <- vector(mode="numeric",length=nrow(fits)-1) #vectors to loop across for BF and cmP
    cmP <- vector(mode="numeric",length=nrow(fits))

  #Loop across sic to compute Bayes factor
  for(i in 1:length(sic)-1){  BF[i] <- exp(sic[i]-sic[i+1]) }
  for(i in 1:length(sic)){    cmP[i] <- (exp(sic[i]-max(sic)))/(sum(exp(sic[1:length(sic)]-max(sic))))    }
    
    
    fits<-  cbind(fits,BF=c(BF,NA),cmP)
  
   return(fits) 
}
