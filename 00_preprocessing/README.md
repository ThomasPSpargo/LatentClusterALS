## 00: Preprocessing scripts

### Polygenic risk scores (PRS)
##### GWAS summary statistic preparation

Publicly available GWAS sumstats for Alzheimer's disease, frontotemporal dementia, Parkinson's disease, and schizophrenia can be obtained and preprocessed by running the `process_sumstatsLCA.sh` script. This script calls to the individual `prep_*_sumstats.R` scripts within the GWAS_sumstats directory and then matches SNP IDs based on chromosomal position to those from a reference panel (see the next section).

Prerequesites for running `process_sumstatsLCA.sh`, including providing the reference panel and the expected directory structure, are described within the script.

Note: Summary statistics used for amyotrophic lateral sclerosis were from a subset of the 2022 GWAS, excluding 'stratum 6' of the GWAS meta-analysis. These are not publically available and were kindly provided by van Rheenen et al. to support our analyses. The only preprocessing of the summary statistics file shared was to rename columns to correspond with the format required to calculate PRS within GenoPredPipe; they were matched to the GenoPredPipe plink-binary reference files using snp (rs)ID. Matching was not performed by chromosomal position since genomic position columns were not included within the file.

For the other GWAS summary statistics, matching by chromosomal position was used in preference to matching by ID to avoid removal of SNPs because of failure to match synonymous rsID marking same genetic variant.


##### PRS calculation and extraction

Polygenic scores were generated via SBayesR using the [GenoPredPipe](https://github.com/opain/GenoPred/tree/master/GenoPredPipe) resource.

The `GenoPredPipe` directory contains the input settings files used with GenoPredPipe.

The script `GenoPredPipe/PRS_extract.R` indicates how PRS scores were extracted from the GenoPredPipe output.

Note that these details are provided for reference but will not be able to run correctly without having available the Project MinE data for which PRS were derived.

### Harmonising Datasets for mplus and subsequent analyses

Data from Project MinE and STRENGTH were preprocessed and harmonised before being entered into Mplus for latent class clustering analysis. The scripts in the `Data_processing` are numbered in order for producing these files. 

An additional file `get_Brainbank_samples.R` is included which can be used to extract the ids for people in the KCL BrainBank.
