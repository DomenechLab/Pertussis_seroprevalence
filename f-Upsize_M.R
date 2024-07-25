#######################################################################################################
# Function to upsize the contact matrix
#######################################################################################################

Upsize_M <- function(M, N, N_u, f_map) {
  #Args: 
  # M: original contact matrix (intensive scale) [square matrix n x n] 
  # N: age-specific population sizes in original population [vector of length n]
  # N_u: age-specific population sizes in upsized population [vector of length n_u > n]
  # f_map: function mapping the augmented age groups to the original age groups [function]
  # Returns: 
  # Upsized matrix Mu [square matrix n_u x n_u]
  
  # Check arguments 
  stopifnot(nrow(M) == ncol(M))
  stopifnot(nrow(M) == length(N))
  stopifnot(length(N_u) > length(N))
  stopifnot(all.equal(sum(N_u), sum(N)))
  n <- length(N)
  n_u <- length(N_u)
  dimnames(M) <- NULL
  
  # Check that matrix M is corrected for reciprocity
  Ni_mat <- matrix(data = NA, nrow = n, ncol = n)
  for(i in 1:n) {
    Ni_mat[i, ] <- N[i] # General term: Ni_mat[i, j] = N[i]
  }
  E <- M * Ni_mat
  stopifnot(isSymmetric(E))
  
  # General term of upsized matrix: M[f(i), f[j]] * N_u[j] / N[f(j)]
  M_u <- matrix(data = NA, nrow = n_u, ncol = n_u)
  for(i in 1:n_u) {
    for(j in 1:n_u) {
      M_u[i, j] <- M[f_map(i), f_map(j)] * N_u[j] / N[f_map(j)]
    }
  }
  
  # Also calculate matrices E and F
  
  # Calculate auxiliary matrices
  Ni_mat <- Nj_mat <- matrix(data = NA, nrow = n_u, ncol = n_u)
  for(i in 1:n_u) {
    Ni_mat[i, ] <- N_u[i] # General term: Ni_mat[i, j] = N_u[i]
    Nj_mat[, i] <- N_u[i] # General term: Nj_mat[i, j] = N_u[j]
  }
  
  E_u <- M_u * Ni_mat
  F_u <- M_u / Nj_mat
  
  stopifnot(isSymmetric(E_u))
  stopifnot(isSymmetric(F_u))
  stopifnot(all.equal(sum(E), sum(E_u)))
  
  out <- list("E_mat" = E_u, "M_mat" = M_u, "F_mat" = F_u)
  
  return(out)
}