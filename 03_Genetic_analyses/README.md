## 03: Genetic analyses

Scripts in this directory are all for analysis of rare and common genetic variation associated with clusters identified by LCA. The analyses are handled by the script `03_genetic_tests.R` which submits the `fishers_compare.R` script once for the rare genetic analysis and performs repeated of `prs_compare.R` with different predictor/data configurations.

The results of each analysis are returned within their own directory. An overview of results from the Fisher's exact tests can be found in the file "Final_summary_tab_\<date\>.tsv"

A figure summarising all PRS analyses can be found in `./PRS_reg/PRS_summaryFig.pdf`.
