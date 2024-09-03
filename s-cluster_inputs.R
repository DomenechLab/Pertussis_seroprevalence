#######################################################################################################
# Script to create text file "s-cluster_inputs.txt" for running simulations on cluster
#######################################################################################################


# Load packages -----------------------------------------------------------
source("s-base_packages.R")
options("useFancyQuotes" = FALSE) # To have standard double quotes
if(file.exists("s-cluster_inputs.txt")) file.remove("s-cluster_inputs.txt")

# List countries ----------------------------------------------------------
count_l <- read.table(file = "_data/list_countries.txt") %>% pluck(1)
rhoV_val <- 0.25
DV_val <- 50
rhoR_val <- 0.66
DR_val <- 66
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