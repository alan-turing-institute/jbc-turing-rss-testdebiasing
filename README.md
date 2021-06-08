## Debiasing targeted testing data

This repository contains the R scripts needed to reproduce the results reported in
the manuscript 'Local prevalence of transmissible SARS-CoV-2 infection : an integrative causal model for debiasing fine-scale targeted testing data'.  

To get started, you will need to install the [renv](https://rstudio.github.io/renv/index.html) R package and run
```
renv::restore()
```
to install the required packages for the scripts.  

### Scripts

`00_download_data.R`  fetches the required data.  
`01_preprocess_data.R`  performs some initial preprocessing of the data for downstream use.  
`02_calculate_infectiousness_estimates` computes the time-dependent probability of being infectiousness conditional on returning a (true) positive PCR test  
`03_prevalence_comparisons.R` calculates prevalence estimates and credible intervals using 1) only Pillar 1+2 data and 2) only REACT data.  
`04a_main_run.R` contains the analysis scripts to generate the results in the main part of the manuscript. 
`04b_main_cut.R` contains the analysis scripts to generate the results for Validation 2. 
`05_supp_run.R` contains the analysis scripts to generate the results for the supplementary material.  
`06*.R` contains the plotting scripts to generate the figures in the main part of the manuscript.  
`07*.R` contains the plotting scripts to generate the figures in the supplementary material.  

`create_SIR_transition_matrices.R` precalculates the hidden Markov model transition matrices in transmats/ needed to fit the stochastic epidemic model.  


### Notes

The script to generate SI Figure 15 uses python. To install the relevant dependencies 
via a conda environment, run `conda env create -f testdebiasing.yml`.