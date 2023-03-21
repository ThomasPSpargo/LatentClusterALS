#Read in files
library(dplyr)

#Read brain bank illumina sample IDS
bb<- tibble(scan("/Users/tom/Downloads/brainbanksamples.txt",what="character"))
names(bb)[1] <- "Illumina_Sample_Barcode"

#Read in the cluster allocations
data <- tibble(read.delim("/Users/tom/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/Writeup/WriteupModels/JNT_JanFinal_savedata_C5/completesample/saved_mplusfit.tsv",sep="\t"))
names(data)[grepl("Illumina_Sample_Barcode",names(data))] <- "Illumina_Sample_Barcode"

#Append cluster assignments
brainb.clust<- data %>%
    dplyr::select(C,Illumina_Sample_Barcode,All_known_IDs) %>%
    right_join(bb,by="Illumina_Sample_Barcode") %>%
  dplyr::filter(!is.na(C))
  
#String split and extract the IDs used in the brain bank
brainb.clust$brainbankIDs<- strsplit(brainb.clust$All_known_IDs,split=",") %>%
  sapply(.,grep,pattern="^A[0-9]+",value=TRUE)

write.table(brainb.clust,file="/Users/tom/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets/brainbank_LCA_clusters.csv",sep=",",row.names = FALSE)
