#######################################################################################################
# Prepare contact matrix for use in serotransmission model
#######################################################################################################

CreateContactMatrix <- function(country_nm, Nvec, source_dat, trim_mat = T, debug = T) {
  
  # Args: 
  # country_nm: name of country [string]
  # Nvec: age-specific population sizes in simulated population [numeric vector]
  # source_dat: source for contact data, either "Mistry" of "Prem" [string]
  # trim_mat: should the SCM be trimmed? If true, the age groups 80-84 yr will be removed [boolean]
  # debug: should the intermediate contact matrices be plotted? [boolean]
  # Returns:
  # 2-list with SCM matrices on density scale (F) and intensive scale (M), per YEAR [2-list matrix]
  
  # Load data
  if(source_dat == "Mistry") {
    path_dat <- "_data/_contact_matrices/Mistry_2021/"
    f_country <- list.files(path = path_dat, pattern = country_nm)
    stopifnot(length(f_country) == 1)
    
    # Load raw data (unit: PER DAY, dimension 85*85, 85 1-yr age groups from age 0 to 84)
    M_mat <- read_csv(file = sprintf("%s/%s", path_dat,f_country), 
                      col_names = F, 
                      show_col_types = F)
    
    if(trim_mat) M_mat <- M_mat[-c(81:85), -c(81:85)] # Trimmed matrix has dimension 80x80
    
    N_tar <- c(sum(Nvec[1:2]), Nvec[-c(1:2)]) # Vector of length 80
    f_map <- function(i) ifelse(i <= 2, 1, i - 1) # Mapping function
  } else if(source_dat == "Prem") {
    
    # Load raw data (unit: PER DAY, dimension 16*16, 16 5-yr age groups from 0-4 to 75-79)
    M_mat <- readRDS("_data/_contact_matrices/Prem_2021/prem_matrices.rds")
    M_mat <- M_mat[[country_nm]]
    
    N_tar <- rep(sum(Nvec) / 16, 16)
    f_map <- function(i) ifelse(i <= 6, 1, ifelse(i >= 77, 16, ceiling((i - 6) / 5) + 1))
  }
  
  # Correct for reciprocity
  colnames(M_mat) <- NULL
  M_mat <- as.matrix(M_mat)
  M_mat_corr <- Project_M(M = M_mat, N_tar = N_tar)
  
  # Upsize contact matrix
  M_u <- Upsize_M(M = M_mat_corr, N = N_tar, N_u = N_vec, f_map = f_map)
  
  # Convert to per year and return
  out <- list("M_mat" = 365 * M_u$M_mat, "F_mat" = 365 * M_u$F_mat)
  
  if(debug) {
    PlotMatrix(M_in = M_mat, plot_title = "Raw M matrix, without correction")
    PlotMatrix(M_in = M_mat_corr, plot_title = "M matrix, corrected for reciprocity")
    PlotMatrix(M_in = M_u$M_mat, plot_title = "M matrix, upsized")
    PlotMatrix(M_in = out$F_mat, plot_title = "F matrix, upsized")
  }
  
  return(out)
}