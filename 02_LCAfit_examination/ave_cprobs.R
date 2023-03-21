suppressPackageStartupMessages({
  library(reshape2) #For the melt function
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

filepath <- commandArgs(trailingOnly = TRUE)

#Obtain directory tree minus model name. This is where output files will be saved
path <- file.path(dirname(filepath),"cprobs/")
if(!dir.exists(path)){dir.create(path,recursive = TRUE)}

#Import data, reading in the cprobs columns and ordering them
data <- data.table::fread(file=filepath, header=TRUE) %>%
  dplyr::select(c("C", paste0("CPROB",1:length(grep("CPROB",colnames(.)))))) 


#Generate probabilities table
prob <- data %>%
  group_by(C) %>%
  summarise(across(starts_with("CPROB"), ~ mean(.x, na.rm = TRUE))) %>%
  round(.,3) %>%
  rename_with(~gsub("CPROB","",.), .cols = everything())


#Save the tabulation in a file
write.table(prob,file=paste0(path,"cprobtable.tsv"),sep="\t",row.names = FALSE,quote=FALSE)
  
  

#Assign group names and then remove first column,containing group levels
prob<- as.matrix(prob)
rownames(prob) <- prob[,1]
prob <- prob[,-1]

#Melt for use in ggplot
mel_cprobs <- reshape2::melt(prob) 

#Convert into factor and assign levels for Var1 and Var2
mel_cprobs$Var1 <- as.factor(mel_cprobs$Var1)
mel_cprobs$Var2 <- as.factor(mel_cprobs$Var2)

#Create ggplot of the fits
plot <- ggplot(data=mel_cprobs,aes(x=Var2,y=Var1,fill=value))+
  geom_tile()+
  geom_text(aes(Var2,Var1, label = formatC(value,format = "f",digits=3)), color = "grey65", size = 4) +
  scale_fill_gradient2(space = "Lab", na.value = scales::muted("blue"), name="Probabilty")+
  theme_bw()+
  xlab("Probability of class")+
  ylab("Assigned class")+
  scale_y_discrete(limits=rev)+
  theme(legend.position = "right")
 
#Save the plot
ggsave(plot=plot, filename = paste0(path,paste0("ave_cprobs_",Sys.Date(),".pdf")),
       units="mm",width=175,height=150)

