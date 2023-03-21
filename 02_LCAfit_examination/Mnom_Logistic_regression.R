library(foreign)
library(nnet)
library(MASS)
library(broom)
library(optparse)

option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file from which a model will be obtained"),
  make_option("--keepcols", action="store", default=NA, type='character',
              help="Indicate columns which will be utilised in the script"),
  make_option("--Predictorcols", action="store", default=NA, type='character',
              help="Indicate columns which will be utilised as model predictors")
)

opt = parse_args(OptionParser(option_list=option_list))

#Extract the columns to keep and to use as predictors
keepcols <- strsplit(opt$keepcols,split=",")[[1]]
predictorcols <- strsplit(opt$Predictorcols,split=",")[[1]]

path <- file.path(dirname(opt$resultsfile),"MnomLogReg/")
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import primary dataset and format
data <- data.table::fread(file=opt$resultsfile, select=keepcols, header=TRUE)
data <-na.omit(data)

tab<- table(data$C)
data$C<- relevel(as.factor(data$C),ref=which(tab==max(tab)))

if("Phenotype" %in% names(data)) data$Phenotype <- factor(data$Phenotype,levels=1:3,labels=c("ALS","PLS","PMA"))

########
#### Run complete case analysis
########

#Run the model
f_rh <- reformulate(predictorcols)
mlr <- multinom(update.formula(f_rh,C~.),data=data)

out<- MASS::stepAIC(mlr,trace = -1) #Run stepwise feature selection


#Summarise the model fit
modsum <- broom::tidy(out)
modsum$beta.l95 <- modsum$estimate-1.96*modsum$std.error
modsum$beta.u95 <- modsum$estimate+1.96*modsum$std.error
modsum$OR <- exp(modsum$estimate)
modsum$l95 <- exp(modsum$beta.l95)
modsum$u95 <- exp(modsum$beta.u95)

modsum$coef95 <- paste0(signif(modsum$estimate,3)," (",signif(modsum$beta.l95,2),", ",signif(modsum$beta.u95,2),")")
modsum$roundedSE <- signif(modsum$std.error,3)
modsum$roundedZ <- signif(modsum$statistic,2)
modsum$roundedP <- signif(modsum$p.value,3)

write.csv(x=modsum, file=paste0(path,"multiLogReg_",Sys.Date(),".csv"))
