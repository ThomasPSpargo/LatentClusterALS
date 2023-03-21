#library(lubridate) #For days to years conversion
library(dplyr) #%>%
library(readxl) #For reading genetic dataset(s)
library(tidyr) #For pivoting

#Specify path containing input files
inpath <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/datasets"

#Specify path out for files
outpath <- "~/OneDrive - King's College London/PhD/PhD project/Latent cluster analysis/final_datasets"
if(!dir.exists(outpath)){dir.create(outpath)}




#-----------------------------------------------#
#             Compute datasets for LCA          #
#-----------------------------------------------#

##################################################
#     Read and preprocess Phenotype datasets     #
#################################################

#Updated Project MinE dataset
NP <- read.csv(file.path(inpath,"New_PM_Pheno2022-04-10.csv"),
               stringsAsFactors = FALSE)

#Read old phenotype data, which includes PC information
PM_phen <- read.delim(file.path(inpath,"recodePM.pheno.2_1_correct.txt"))

#Subset NP data to be only those who were in the old dataset, and thus QC protocol completed
old_sample <- which(NP$Illumina_Sample_Barcode %in% PM_phen$Illumina_Sample_Barcode) #Find matching barcodes
NP.wgs <- NP[old_sample,]

  ################################################
  #    Read and preprocess  Genetic datasets     #
  ################################################
  
##Read in gene-panel genetic data and merge into one data.frame
  gene_vars <- read_excel(file.path(inpath,"results_puja_new_list_of_genes.xlsm"), 
                         na="NA", sheet = "variants")
  gene_vars_c9 <- read_excel(file.path(inpath,"results_puja_new_list_of_genes.xlsm"), 
                            na="NA", sheet = "c9 status")
  gene_vars_atxn2 <- read_excel(file.path(inpath,"results_puja_new_list_of_genes.xlsm"), 
                               na="NA", sheet = "atxn2 status")
  
  colnames(gene_vars_atxn2)[which(colnames(gene_vars_atxn2)=="atxn2 status")] <- "atxn2_status"   #Tweak coding to remove space
  #Rename ID column for comparability - for merging
  colnames(gene_vars)[1] <- colnames(gene_vars_c9)[1] <- colnames(gene_vars_atxn2)[1] <- colnames(NP)[1]
  
  
##READ in PRS data
  gene_vars_prs <- read.csv(file.path(inpath,"SBayesR_PRS.csv"),
             stringsAsFactors = FALSE)
  
  
  #####
  ### Collapse gene data across moderate and high impact variants
  #####
  
  #Extract relevant gene names and then the 'mod' labels, to avoid doubling up genes; sub out additional string portions
  col <- colnames(gene_vars)[3:ncol(gene_vars)] %>%
    grep("mod", ., value = T) %>%
    gsub("(_|.)mod","", .)
  
  #Set function to identify total numbers of variants each person had in each gene, aggregating across high and mod impact columns
    #Col is a vector containing patterns indicating the columns of interest  (i.e. gene name)
    #allCols is the full dataframe to search across,
  aggregateGenes <- function(col,allCols){
    
    #Identify columns containing the pattern to match
    gen_index<- grep(col,colnames(allCols)[1:ncol(allCols)]) 
    
    #Generate which aggregates across the two elements selected, if all is NA, then set NA, otherwise rowsum
    out<- allCols %>%
      select(all_of(gen_index)) %>%
      mutate(joint=case_when(if_any(,~!is.na(.)) ~ rowSums(across(everything()),na.rm = TRUE),
                        if_all(,~is.na(.)) ~ NA_real_
        )) %>%
      select(3)
    
    #Return as simple vector
    out <- unname(as.vector(out))
    
    return(out)
  }
    
  #Sapply across all gene names, and cbind to overwrite the original gene_vars
  genes_total<- sapply(col,aggregateGenes,allCols=gene_vars[-(1:2)])
  gene_vars <- cbind(gene_vars[,1:2],genes_total)
  
  ######################################
  # Combine pheno and genetic datasets #
  ######################################
  
  #Combine all datasets
  PMFULLdat <- NP.wgs %>%
    full_join(gene_vars,by = colnames(NP.wgs)[1]) %>%       #Binary vars
    full_join(gene_vars_c9,by = colnames(NP.wgs)[1]) %>%    #C9 status
    full_join(gene_vars_atxn2,by = colnames(NP.wgs)[1]) %>%   #Atxn2 status
    left_join(gene_vars_prs,by = colnames(NP.wgs)[1]) %>%     #PRS scores - left join bring across only people from the original data freeze 
    left_join(PM_phen[c("Illumina_Sample_Barcode",grep("^PC[0-9]+",colnames(PM_phen),value=TRUE))],by="Illumina_Sample_Barcode") #Genetic principal components
  
  
  ######################################
  #         QC combined dataset        #
  ######################################
  
  #####
  ###  Filter data from cases and controls in combined dataset
  #####
  
  #Save into new object to avoid overwriting
  ##FILTERING of cases
  #Drop controls and people with 'other' diagnosis
  #Keep only unrelated individuals
  #Retain only people passing QC
  
  #Extractthe control cohort, dropping any cols which are NA for all people. Intent is to compare prs to a control reference rather than just within clusters.
  PM_control <- PMFULLdat %>%
    filter(Phenotype=="Control" &
             relatedness=="unrelated_individual" &
             failed_QC=="N") %>%
    .[,colSums(is.na(.))<nrow(.)]
  
  # # #Sample size/Exclusion counts:
  # PMFULLdat %>%
  #   filter(Phenotype=="Control") %>%
  #   group_by(relatedness) %>%
  #   summarise(n())
  #
  # PMFULLdat %>%
  #   filter(Phenotype=="Control" &
  #            relatedness=="unrelated_individual") %>%
  #   group_by(failed_QC) %>%
  #   summarise(n())
  #
  # PMFULLdat %>%
  #   filter(Phenotype=="Control" &
  #            relatedness=="unrelated_individual" &
  #            failed_QC=="N") %>%
  #   summarise(n())
  
  #Write the control data
  write.table(PM_control,file=file.path(outpath,"PM_CONTROLS_data.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)
  
  ## Extract the case data
  PMFULLfilt <- PMFULLdat %>%
    filter(!Phenotype %in% c("Control","other","FTD", NA_character_) &
             relatedness=="unrelated_individual" &
             failed_QC=="N")
  
  # # #Sample size/Exclusion counts:
  # PMFULLdat %>%
  #   filter(!Phenotype %in% c("Control","other","FTD", NA_character_)) %>%
  #   group_by(relatedness) %>%
  #   summarise(n())
  #
  # PMFULLdat %>%
  #   filter(!Phenotype %in% c("Control","other","FTD", NA_character_) &
  #            relatedness=="unrelated_individual") %>%
  #   group_by(failed_QC) %>%
  #   summarise(n())
  #
  # PMFULLdat %>%
  #   filter(!Phenotype %in% c("Control","other","FTD", NA_character_) &
  #            relatedness=="unrelated_individual" &
  #            failed_QC=="N") %>%
  #   summarise(n())
  
  ######
  ### Mutate columns in the case data
  ######
  
  ##Phenotype: ALS variants to ALS specific disease; Retain FTD as separate - clump ALS/FTD because it wont be a clean diagnosis (other 'pure' ALS diagnoses will have non-recorded FTD features)
  ## atxn2_status: make binary, any expansion considered variant
  ## c9_status: make binary. Inconsistent represents technical failures, intermediate (<30 repeats, treat as no variant)
  PMFULLfilt <- PMFULLfilt %>%
    mutate(Phenotype = case_when(grepl("ALS",Phenotype) ~ "ALS",
                                Phenotype == "PBP" ~ "ALS",
                                TRUE ~ Phenotype),
           atxn2_status = case_when(atxn2_status=="normal" ~ 0,
                                  atxn2_status %in% c("intermediate","expanded") ~ 1,
                                  TRUE ~ NA_real_),
           c9_status = case_when(c9_status %in% c("normal","intermedidate") ~ 0,
                               c9_status=="inconsistent" ~ NA_real_,
                               c9_status== "expanded" ~ 1,
                               TRUE ~ NA_real_),
           new_num_id = as.numeric(paste0(999,row_number())) #Assign unique numerical id number - important for matching mplus output with the raw file
    )
  
  
  ######
  ###   Recode survival impossible values in survival and disease duration variables
  ######
  
  #Identify people with delay of 0 or less 
  del_0 <- which(PMFULLfilt$Diagnostic_delay_years<=0)
  PMFULLfilt$Diagnostic_delay_years[del_0] <- NA #Delay of 0 is not possible, code as missing
  
  #Some people have a reported 0 days of survival. This is impossible.
  #Recode the people with survival of 0 days as having at least 1 day of survival but are censored
  surv_0 <- which(PMFULLfilt$Time_to_death_or_last_followup_years<=0)
  PMFULLfilt$Time_to_death_or_last_followup_years[surv_0] <- 1/365.25 #Assign 1 day, on years scale
  PMFULLfilt$survival_status_bin[surv_0] <- 0                         #The person is censored
  
  #Add-in records of people with Diag Delay info but no survival, then flag them as censored
  missing_surv<- which(!is.na(PMFULLfilt$Diagnostic_delay_years) & is.na(PMFULLfilt$Time_to_death_or_last_followup_years))
  PMFULLfilt$Time_to_death_or_last_followup_years[missing_surv] <- PMFULLfilt$Diagnostic_delay_years[missing_surv]
  PMFULLfilt$survival_status_bin[missing_surv] <- 0
  
  
  ####Next, check for mismatched missingness in SRV and SRVBIN,
  #code any unpaired records as NA if SRV is missing or censored if SRVBIN is missing.
  PMFULLfilt$survival_status_bin[which(!is.na(PMFULLfilt$Time_to_death_or_last_followup_years) & is.na(PMFULLfilt$survival_status_bin))] <- 0
  PMFULLfilt$survival_status_bin[which(is.na(PMFULLfilt$Time_to_death_or_last_followup_years) & !is.na(PMFULLfilt$survival_status_bin))] <- NA
  
  
  ################################################
  ###  Generate genetic aggregation columns    ###
  ################################################
    
  #String labels for genes by functional groups
  auto.prot   <- c("UBQLN2","CHMP2B","ERBB4","SIGMAR1","CHCHD10","DAO","OPTN","SCFD1","SQSTM1","TBK1","UNC13A","VAPB","VCP","VEGFA")
  RNA.func    <-  c("ANG","atxn2_status","FUS","HNRNPA1","MATR3","SETX","TAF15","TARDBP")
  Cyto.transp <- c("ALS2","ANXA11","C21ORF2","FIG4","NEFH","NEK1","ATXN1","DCTN1","MOBP","PFN1","SPG11","TUBA4A")
  
  #Create also an indices of every relevant gene column across pathways, and then with C9 and SOD1
  allg <- c(auto.prot,RNA.func,Cyto.transp)
  allgSOD1C9 <- c(allg,"SOD1","c9_status")
  
  #Generate and binary-code aggregate columns
  PMFULLfilt <- PMFULLfilt %>%
    mutate(across(all_of(allg), ~ case_when(. > 1 ~ 1, TRUE ~ .)), ### FIRSTLY, recode any 'multiple' variants in one gene to 1
           AutoProt=case_when(if_any(all_of(auto.prot),~!is.na(.)) ~ rowSums(across(all_of(auto.prot)),na.rm = TRUE),
                                  if_all(all_of(auto.prot),~is.na(.)) ~ NA_real_
           ),
           AutoProt=if_else(AutoProt>=1,1,AutoProt),
           RNAFunc=case_when(if_any(all_of(RNA.func),~!is.na(.)) ~ rowSums(across(all_of(RNA.func)),na.rm = TRUE),
                              if_all(all_of(RNA.func),~is.na(.)) ~ NA_real_
           ),
           RNAFunc=if_else(RNAFunc>=1,1,RNAFunc),
           CytoTransp=case_when(if_any(all_of(Cyto.transp),~!is.na(.)) ~ rowSums(across(all_of(Cyto.transp)),na.rm = TRUE),
                             if_all(all_of(Cyto.transp),~is.na(.)) ~ NA_real_
           ),
           CytoTransp=if_else(CytoTransp>=1,1,CytoTransp),
           BurdenG=case_when(if_any(all_of(allg),~!is.na(.)) ~ rowSums(across(all_of(allg)),na.rm = TRUE),
                                if_all(all_of(allg),~is.na(.)) ~ NA_real_
           ),
           BurdenG=if_else(BurdenG>=1,1,BurdenG),
           )
    
# ######
# ### Write genetic summary to file - NOT RUN
# ######
# geneSelect <- PMFULLfilt[c("new_num_id",
#                            "c9_status",
#                            "SOD1",
#                            "AutoProt",
#                            "RNAFunc",
#                            "CytoTransp",
#                            "BurdenG",
#                            allg,
#                            grep("(SBayesR|^PC[0-9]+)",colnames(PMFULLfilt),value=TRUE))]
#                               
# 
# ## Write genetic variables file
# write.table(geneSelect,file=file.path(outpath,"PM_filtered_all_genetics.csv"),row.names=FALSE, col.names = TRUE,sep=",",quote = FALSE)

  ######
  ### Write index file for identifying samples 
  ######
  
  PMFULLfilt[c("new_num_id",
               "Illumina_Sample_Barcode",
               "Patient_ID",
               "All_known_IDs",
               "origin")] %>%
    write.table(.,file=file.path(outpath,"PM_identities.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)
  
  
  ####Transform chr columns into numeric, to be parsed by mplus
  
  PMFULLfilt <- PMFULLfilt %>%
    mutate(
      Sex_at_birth = case_when(Sex_at_birth == "male" ~ 1,
                               Sex_at_birth == "female" ~ 2
      ),
      Phenotype = case_when(Phenotype == "ALS" ~ 1,
                        Phenotype == "PLS" ~ 2,
                        Phenotype == "PMA" ~ 3,
      ),
      Site_of_Onset = case_when(Site_of_Onset == "spinal" ~ 1,
                         Site_of_Onset == "bulbar" ~ 2,
                         Site_of_Onset == "generalized" ~ 1,
                         Site_of_Onset == "thoracic/respiratory" ~ 1,
                         Site_of_Onset == "other" ~ 1, #Other category present within STR data
                         TRUE ~  NA_real_
      )
  ) 
  
  
  #Write file containing all fields from the dataset
  write.table(PMFULLfilt,file=file.path(outpath,"PM_LCA_Allfields.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)
  
  
  
  ######
  ### Generate a PM-dataset per-country standardisation panel for diagnostic delay
  ######
  #Generate a per-origin standardisation panel for country-wise diagnostic delay
  delay_panel<- PMFULLfilt %>%
    filter(!is.na(Diagnostic_delay_years)) %>%
    group_by(origin) %>%
    summarise(N=n(),
              mean=mean(Diagnostic_delay_years),
              sd=sd(Diagnostic_delay_years)
    ) %>%
    mutate(ID=origin) %>%
    rename(iso2c=origin)
  
  write.table(delay_panel,file=file.path(outpath,"PM_diagDelay_country_standardised.tsv"),row.names=FALSE, col.names = TRUE,sep="\t",quote = FALSE)