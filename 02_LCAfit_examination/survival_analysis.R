suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(optparse)
})

option_list = list(
  make_option("--resultsfile", action="store", default=NA, type='character',
              help="Specify the source file from which a model will be obtained"),
  make_option("--coeftableFunpath", action="store", default=NA, type='character',
              help="File path to function for parsing Cox coefficient's table")
)

opt = parse_args(OptionParser(option_list=option_list))

#Source custom function for pretty-generation of a coefficients table
source(opt$coeftableFunpath)

#Obtain directory tree minus model name with appended file. This is where output files will be saved
path <- file.path(dirname(opt$resultsfile),"SurvivalAnalysis")
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import dataset
data <- tibble(data.table::fread(file=opt$resultsfile,header=TRUE))

#Palette for plot
cbbPalette <- c("gold", "violetred", "#56B4E9", "#000000", "limegreen", "#D55E00", "#CC79A7", "#F0E442")

######
### Univariate analysis
######

#Fit survival curve
KMfit<- survfit(Surv(Time_to_death_or_last_followup_years,survival_status_bin) ~ C,data=data)
names(KMfit$strata)<- gsub("C=","", names(KMfit$strata))

#Build KM plot
kmplot <- ggsurvplot(KMfit,data=data,
                    pval = F,pval.method=F, conf.int = T, risk.table = F,
                    xlab= "Time (years)",
                    legend="right",
                    palette = cbbPalette,
                    legend.title="Class"
)

#Run pairwise logrank tests
diff <- pairwise_survdiff(Surv(Time_to_death_or_last_followup_years,survival_status_bin) ~ C,data=data, p.adjust.method = 'fdr')

#Save the figure
ggsave(plot=kmplot[[1]], 
       filename = file.path(path,"KM_plot.pdf"),
       units="mm",width=200,height=150,dpi = 300)

#Save the pairwise survdiff results
data.table::fwrite(data.frame(diff$p.value),
                   file=file.path(path,"survdiff.csv"),
                   sep=",",row.names = TRUE,col.names = TRUE)

#Assign reader-friendly names for the predictor columns
predictorcols <- c(Class="C",
                   `Diagnostic Delay`="CTYDEL",
                   `Age of onset`="AGEONS_nml",
                   `Site of onset`="Site_of_Onset",
                   Sex="Sex_at_birth",
                   Diagnosis="Phenotype",
                   `Site of onset`="Site_of_Onset")

#Convert appropriate vars to factor
dcox<- data %>%
  mutate(across(c(C,origin),~as.factor(.)),
         Site_of_Onset=factor(Site_of_Onset,levels=1:2,labels=c('Other','Bulbar')),
         Sex_at_birth=factor(Sex_at_birth,levels=1:2,labels=c('Male','Female')),
         Phenotype=factor(Phenotype,levels=1:3,labels=c("ALS","PLS","PMA"))
         ) %>%
  rename(all_of(predictorcols))

######
### Cox (multivariate) analysis
######
f_rh <- reformulate(paste0("`", names(predictorcols), "`", collapse = "+"))
Coxfit  <- coxph(update.formula(f_rh,Surv(Time_to_death_or_last_followup_years,survival_status_bin)~.), data=dcox)

#Apply custom function to format cox coefficients table
coefTab<- coxCoefTable(Coxfit,dcox,cutpoints=NULL,truncatecutpoints=NULL)

#Save the forest plot
ggsave(plot=coefTab$figure, filename = file.path(path,"COXforest.pdf"),
       units="mm",width=200,height=150,dpi = 300)

#Save the tabular version of the coefficients
data.table::fwrite(coefTab$tabulation,
                   file=file.path(path,paste0("COXcoefs_",Sys.Date(),".csv")),
                   sep=",",row.names = FALSE,col.names = TRUE
)
