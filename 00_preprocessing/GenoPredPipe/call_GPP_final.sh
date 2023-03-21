#!/bin/bash
#To be called with GenoPredPipe as the working directory

#Run the main pipeline
snakemake --profile slurm --use-conda run_create_reports --latency-wait 20

#run sbayesr for ALS
snakemake --profile slurm --use-conda resources/data/target_checks/ALS_ProjectMinE/target_prs_sbayesr_EUR_ALScc.done --latency-wait 20

#run sbayesr for schizophrenia
snakemake --profile slurm --use-conda resources/data/target_checks/ALS_ProjectMinE/target_prs_sbayesr_EUR_SZposMatch.done --latency-wait 20

#run sbayesr for Parkinson's Disease
snakemake --profile slurm --use-conda resources/data/target_checks/ALS_ProjectMinE/target_prs_sbayesr_EUR_PDposMatch.done --latency-wait 20

#run sbayesr for Alzheimer's Disease
snakemake --profile slurm --use-conda resources/data/target_checks/ALS_ProjectMinE/target_prs_sbayesr_EUR_AZposMatch.done --latency-wait 20

#run sbayesr for frontotemporal dementia
snakemake --profile slurm --use-conda resources/data/target_checks/ALS_ProjectMinE/target_prs_sbayesr_EUR_FTDposMatch.done --latency-wait 20
