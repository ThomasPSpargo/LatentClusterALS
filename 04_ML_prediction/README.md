## 04: Machine learning prediction

Machine learning classification algorithms (random forest and xgboost) are trained using scripts within this directory.

### Training ML algorithms and general interpretation

Running `04a_ML_predictions.sh` generates the datasets (containing different cohort subsamples and features to train upon) to be supplied to the ML algorithms via `DefineMLdatasets.R`. Once datasets are generated, it then submits slurm scheduler jobs to run R scripts that train and optimise random forest (`random_forest_job.R`) and xgboost (`xgboost_stagedTuning.R`) algorithms for each generated dataset. Models trained for each dataset are returned within their own respective directory. Various summary info about the algorithm performance is obtained via the `ML_extract.R` script that is called automatically when the best-fit model is identified.


### Generating SHAP values

SHapley Additive exPlanations (SHAP) are a model agnostic approach to examining how model features contribute to predictions made by machine-learning algorithms.

`04b_SHAP_predictions_jobs.sh` submits slurm jobs which run the R script `SHAP_job.R`. This script computes SHAP values via [kernelshap](https://github.com/ModelOriented/kernelshap) (v0.3.3) for a given ML models. Models to analyse should be indicated within a plain-text listing file paths to the .Rds files output from jobs generated when running `04a_ML_predictions.sh` (see script for details). The SHAP values returned will be saved within an .Rds output file that can be loaded into R and extracted accordingly. My extraction script is part of `Additional_plots_misc.Rmd`.
