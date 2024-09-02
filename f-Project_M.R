#######################################################################################################
# Project contact matrix matrix to target population structure
#######################################################################################################

Project_M <- function(M, N_tar) {
  # Args: 
  # M: contact matrix in original population [square matrix n x n] 
  # N_tar: vector of age-specific population sizes in target population [vector of length n]
  # Returns: 
  # Matrix nxn, corrected for reciprocity in the target population
  
  #Check arguments 
  stopifnot(nrow(M) == ncol(M))
  stopifnot(nrow(M) == length(N_tar))
  n <- nrow(M)
  
  # Create matrix Nij = Ni
  Ni_mat <- matrix(data = NA, nrow = n, ncol = n)
  for(i in 1:nrow(Ni_mat)) {
    Ni_mat[i, ] <- N_tar[i]
  }
  
  # General term of corrected matrix: (Mij * Ni + Mji * Nj) / (2 * Ni)  
  # Method M1 (pairwise correction) in Ref 1
  E <- M * Ni_mat
  M_out <- (E + t(E)) / (2 * Ni_mat)
  
  return(M_out)
}