#NOTES
#SVALUES instruction here
  #http://www.statmodel.com/discussion/messages/13/3657.html?1471304009

#Function to run mplus and obtain savedata for the model indicated - to be used on the best fit model of each class
suppressMessages(library(MplusAutomation))
library(data.table)
library(optparse)

option_list = list(
  make_option("--Mplusoutfile", action="store", default=NA, type='character',
              help="Specify the mplus model from which to save data"),
  make_option("--PreMplusFields", action="store", default=NA, type='character',
              help="Specify a dataset to merge with the MPLUS output, taking minimal columns from MPlus (Cprobs and class assignment)"),
  make_option("--PreMplusIDname", action="store", default="ID", type='character',
              help="Specify the column name in the pre mplus dataset upon which to match, defaults to 'ID'"),
  make_option("--MplusKeepcolumns", action="store", default="*", type='character',
              help="Pattern passed to grepl, indicates which columns from the Mplus results should be retained. Defaults to '*', which returns all columns. Use | to indicate OR for multiple patterns"),
  make_option("--helperFunDir", action="store", default=".", type='character',
              help="Directory containing custom helper functions leveraged in the script")
)

opt = parse_args(OptionParser(option_list=option_list))

#Source helper functions
invisible(lapply(list.files(opt$helperFunDir,full.names = TRUE,pattern=".R"),source))

#Read in the model and prepare path and prepare directory structure
mod <- readModels(opt$Mplusoutfile)
out_dir <- file.path(dirname(opt$Mplusoutfile),paste0("savedata_C",mod$summaries$NLatentClasses))
if(!dir.exists(out_dir)){dir.create(out_dir,recursive = TRUE)}

message("outdir: ", out_dir)
#Make the file output location
out<- file.path(out_dir,gsub(".out","savedata.inp",mod$summaries$Filename))

#Create a reusable file name for the data output, substring this if too long
dataout <- file.path(out_dir,gsub("_c.*",".dat",basename(out)))

#Convert the model into an mplusAutomation mplusObject with custom function and write the datafile to the output directory
prep <- prepObject(mod,out,dataout)
obj <- prep$obj
dataout <- prep$dataout

obj$TITLE <- paste0(obj$TITLE,"_savedata")   #Adjust object title

obj$SAVEDATA <- paste0("FILE IS ", obj$TITLE,".dat;\n",
                       "FORMAT IS FREE;\n",
                       "SAVE = CPROB;\n",
                       "ESTIMATES = ",obj$TITLE,"_params.txt;") #Set option to save data  

obj$OUTPUT <- paste0(obj$OUTPUT,"\n TECH1;\n",
                     "ENTROPY;\n",
                     "SVALUES;")

#Adjust starts to 0
obj$ANALYSIS <- gsub(paste0("STARTS = ",mod$input$analysis$starts, ";"), 
                     paste0("STARTS =  0;"),
                     obj$ANALYSIS)

if(mod$input$analysis$starts!=0){ #If the model does not already have 0 starts set (indicating the fit has been predefined)
  
  #Identify row of output containing the seed values
  seedpos<- grep("Final stage loglikelihood values at local maxima, seeds, and initial stage start numbers:", mod$output)+2
  seed <- as.numeric(strsplit(mod$output[seedpos]," ")[[1]]) #String split the best seed row, and coerce to numeric - blank space will become NA
  seed <- seed[!is.na(seed)][2] #Subset to only non-na values, and take the 2nd record, the seed value
  
  obj$ANALYSIS <- paste0(obj$ANALYSIS,"\nOPTSEED = ", seed, ";")   #set seed number
}

#Run save data model
mod2 <- mplusModeler(obj, modelout=out,
                     dataout = dataout,
                     writeData='never', #This is the only allowed option; set to suppress a warning message
                     hashfilename = TRUE, #FALSE in order to avoid 90 character file limits
                     run=1L
)

#Write out fit for sanity check
message("Loglikelihood for model input: ",mod$summaries$LL,". LL from savedata output: ",mod2$results$summaries$LL)

#Extract the saved dataset
data <- mod2$results$savedata

#Optionally merge with additional fields not used by mplus
if(!is.na(opt$PreMplusFields)){
  fullDF<- fread(opt$PreMplusFields)
  
  colnames(data)[colnames(data)=="ID"] <- opt$PreMplusIDname
  keepcols<- c(opt$PreMplusIDname,grep(opt$MplusKeepcolumns,colnames(data),value=TRUE))
  
  message("Joining datasets and Keeping the mplus columns: ",paste(keepcols,collapse= ", "))
  data <- left_join(data[,keepcols],fullDF,by=opt$PreMplusIDname)
}

#Write out
finalfile <- file.path(out_dir, "saved_mplusfit.tsv")
data.table::fwrite(data,file=finalfile,sep="\t")
message(finalfile," file written")

#Print output file name to console
cat(finalfile)
