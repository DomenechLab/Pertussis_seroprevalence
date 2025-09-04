#######################################################################################################
# Script to create text file "s-cluster_inputs.txt" for running simulations on cluster
#######################################################################################################


# Load packages -----------------------------------------------------------
source("s-base_packages.R")
options("useFancyQuotes" = FALSE) # To have standard double quotes
if(file.exists("s-cluster_inputs.txt")) file.remove("s-cluster_inputs.txt")
if(file.exists("s-cluster_inputs2.txt")) file.remove("s-cluster_inputs2.txt")

# List countries for simulations at equilibrium ----------------------------------------------------------
count_l <- ct_list
rhoV_val <- rhoR_val <- 2
DV_val <- DR_val <- 50
v_cov_val <- 0.9

for(s in count_l) {
  cat("sbatch --export=COUNTRY=", dQuote(s), 
      ",RHO_V=", dQuote(rhoV_val), 
      ",D_V=", dQuote(DV_val), 
      ",RHO_R=", dQuote(rhoR_val), 
      ",D_R=", dQuote(DR_val), 
      ",V_COV=", dQuote(v_cov_val), 
      " m-run_simulations.sh \n",
      sep = "", 
      append = T, 
      file = "s-cluster_inputs.txt")  
}

# List countries for simulations and empirical comparison ----------------------------------------------------------
count_l <- read.table(file = "_data/list_countries2.txt") %>% pluck(1)
rhoV_val <- rhoR_val <- 5
DV_val <- DR_val <- 90

for(s in count_l) {
  cat("sbatch --export=COUNTRY=", dQuote(s), 
      ",RHO_V=", dQuote(rhoV_val), 
      ",D_V=", dQuote(DV_val), 
      ",RHO_R=", dQuote(rhoR_val), 
      ",D_R=", dQuote(DR_val), 
      " m-run_simulations_empirical_comparison.sh \n",
      sep = "", 
      append = T, 
      file = "s-cluster_inputs2.txt")  
}