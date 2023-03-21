## 01: Fitting LCA

### NOTES 
- Please update file paths in the scripts `01a_create_inp_files_plus_bash_scripts.sh`, `01b_submit_mplusjobs.sh`, `iterateModels.R`, `template_mplusjob.tpl`, `graph_IC.R`.

- LCA is a maximum-likelihood approach and Mplus identifies the best fit for the provided data and a given number of classes by specifying random starting values and iteration until the best fit is obtained. The best fit model for a given number of classes in the provided data is indicated by replication of the smallest absolute loglikelihood value from multiple random starts.

- Best-fitting latent class models are identified using the `iterateModels.R` script which uses the MplusAutomation R package and custom helper functions provided in this repository to run an initial Mplus input file and then subsequent Mplus input files with adjusted options (e.g., the `STARTS` and `STSCALE` commands) that are generated when running the script until the best-fiting model is identified for a given number of classes.

### Preparing for latent class analysis (LCA)

Run `bash 01a_create_inp_files_plus_bash_scripts.sh` to generate initial Mplus input files from the plain text file `mplus_inp_template` and corresponding bash scripts to run the analyis from the `template_mplusjob.tpl` template. Bash scripts are output in the same directory as each set of Mplus input files.

The format of `mplus_inp_template` is parsed by the [MplusAutomation R package](https://cran.r-project.org/web/packages/MplusAutomation/vignettes/vignette.html) (version 1.1.0) and converted into a series of initial mplus input files for clustering with different dataset configurations and specifies to attempt models of 1 to 9 classes.

The `template_mplusjob.tpl` template specifies to a bash script which indicates to run `iterateModels.R` across all mplus input files in a given directory (henceforth 'all models') and then, after the best fit is found for all models, to run the script `graph_IC.R` compares several information criterion metrics across the models fitted for each number of classes.

### Fitting LCA

Once the analysis files have been generated, each bash script can be submitted to a slurm scheduler using: `bash 01b_submit_mplusjobs.sh`.

This submission script checks for expected outputs from each potential job and will only submit the job file if an output is missing. (e.g., jobs will not be submitted if a `bestout` subdirectory exists (see below) and contains Mplus `.out` files for all models.)

### Outputs

For each set of models fitted, outputs are returned across several directories.

Each initial Mplus input (.inp) file will be accompanied by a corresponding .out file that indicates the model fitted using the specifications from `mplus_inp_template`. Subsequent .inp and .out files generated when running `iterateModels.R` are returned in the subdirectory `allfits`. The .inp and .out files which specify the best fitting model for a given number of classes are copied from `allfits` into the directory `bestout`. The `bestout` directory contains the best-fitting model identified for each number of classes.

The `bestout` directory also contains files which summarise the best-fitting models for each number of classes to give an initial overview of the solution for each number of classes. The individual .out files should also be checked to ensure that a given solution has converged correctly and the model is identified.
