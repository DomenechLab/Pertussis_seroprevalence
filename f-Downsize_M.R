#######################################################################################################
# Function to downsize contact matrix
#######################################################################################################

Downsize_M <- function(M, N, N_d, f_map) {
  #Args: 
  # M: original contact matrix (intensive scale) [square matrix n x n] 
  # N: age-specific population sizes in original population [vector of length n]
  # N_d: age-specific population sizes in downsized population [vector of length n_d < n]
  # f_map: function mapping the downsized age groups to the original age groups [function]
  # Returns: 
  # Upsized matrix M_d [square matrix n_d x n_d]
  
  # Check arguments 
  stopifnot(nrow(M) == ncol(M))
  stopifnot(nrow(M) == length(N))
  stopifnot(length(N_d) < length(N))
  stopifnot(sum(N_d) == sum(N))
  n <- length(N)
  n_d <- length(N_d)
  
  dimnames(M) <- NULL
  
  # Calculate matrix on extensive scale
  Ni_mat <- matrix(data = NA, nrow = n, ncol = n)
  for(i in 1:n) {
    Ni_mat[i, ] <- N[i] # General term: Ni_mat[i, j] = N[i]
  }
  E <- M * Ni_mat
  #browser()
  stopifnot(isSymmetric(E))
  
  # General term of upsized matrix on extensive scale: E_d[i, j] = sum_{i \in f(i), f \in f(j)} Eij
  E_d <- matrix(data = NA, nrow = n_d, ncol = n_d)
  for(i in 1:n_d) {
    for(j in 1:n_d) {
      E_d[i, j] <- sum(E[f_map(i), f_map(j)])
    }
  }
  
  # Calculate matrices M and F
  
  # Calculate auxiliary matrices
  Ni_mat <- Nj_mat <- matrix(data = NA, nrow = n_d, ncol = n_d)
  for(i in 1:n_d) {
    Ni_mat[i, ] <- N_d[i] # General term: Ni_mat[i, j] = N_d[i]
    Nj_mat[, i] <- N_d[i] # General term: Nj_mat[i, j] = N_d[j]
  }
  
  M_d <- E_d / Ni_mat
  F_d <- M_d / Nj_mat
  
  stopifnot(isSymmetric(E_d))
  stopifnot(isSymmetric(F_d))
  stopifnot(sum(E) == sum(E_d))
  
  out <- list("E_mat" = E_d, "M_mat" = M_d, "F_mat" = F_d)
  
  return(out)
}