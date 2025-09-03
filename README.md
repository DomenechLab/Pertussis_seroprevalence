# Codes for pertussis serotransmission model

This folder contains R and C codes to reproduce the results described in the study, "Natural immune boosting biases pertussis infection estimates in seroprevalence studies."

## Required packages 

The model was developed using the package "pomp" (version 5.10) in R version 4.4.1. Other required packages are listed in the script "s-base_packages.R"

## Main scripts

- The model is created via a call to the function "f-CreateSerotransModel.R", which internally calls the C file "c-serotransmission_model.c"
- An example model simulation can be run using the script "m-check_model.R"
- The clustering of next-generation matrices (NGMs) is performed via the script "m-analyze_SCMs.R" 
- The parameter estimation in six European countries is done via the script "m-run_simulations_empirical_comparison.R"
- The model predictions in twelve countries are run using the script "m-run_simulations.R"
- The predictions for the homogeneous model with an alternative structure and parametrization of immune boosting are run using the script "m-run_simulations_alt_model_homogeneous.R" 
- All other scripts are either auxiliary functions (prefix f-), support scripts (prefix s-), or scripts to make the figures

## Data sources 

- Social contact matrix (SCM) and demographic data:
  - [Mistry et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33436609/)
  - Associated online repository: https://github.com/mobs-lab/mixing-patterns 
- Seroprevalence data:
  - [Pebody et al., 2005](https://pubmed.ncbi.nlm.nih.gov/15724723/) (Appendix table)
  - [Wehlin et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34120372/) (Table 3)
