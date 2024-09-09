#######################################################################################################
# Compute R0 for the age-structured model at the disease-free equilibrium without vaccination
#######################################################################################################

compute_R0 <- function(q, alphaV = 0, epsilonV = 0, theta, v1 = 0, epsA = 0, N, delta, Cmat, gamma, type = "SIR") {
  # Args:
  # q: susceptibility factor
  # alphaV: waning rate, per year
  # epsilonV: leakiness of vaccine-derived immunity
  # theta: relative infectiousness of repeat infections
  # v1: vaccine coverage for primary immunization
  # epsA: primary vaccine failure
  # N: vector of population sizes
  # delta: vector of aging rates, per year
  # Cmat: contact rate matrix (per year)
  # gamma: recovery rate
  # type: "SIR", basic model with no repeat infection, or "S2I2R" with repeat infections
  # Returns: next generation matrix
  
  stopifnot(type %in% c("SIR", "S2I2R"))
  kron <- function(i, j) return(ifelse(i == j, 1, 0)) # Kronecker delta function
  N <- unname(N)
  nA <- length(N)
  F_mat <- V_mat <- matrix(0, nrow = 2 *  nA, ncol = 2 * nA) # Create matrices
  
  p_age <- delta / (delta + alphaV) # Probability of aging before immunity wanes 
  
  # First compute equilibria at disease-free equilibrium 
  S1 <- S2 <- V <- numeric(nA)
  S1 <- N # Disease-free equilibrium without vaccination
  theta <- 0

  # Compute the matrices required to calculate the next generation matrix
  F11 <- F12 <- F21 <- F22 <- matrix(0, nrow = nA, ncol = nA)
  V11 <- V12 <- V21 <- V22 <- matrix(0, nrow = nA, ncol = nA)
  
  for(i in 1:nA) {
    for(j in 1:nA) {
      F11[i, j] <- q[i] * S1[i] * Cmat[i, j] / N[j]
      F21[i, j] <- q[i] * (S2[i] + epsilonV * V[i]) * Cmat[i, j] / N[j]
      V11[i, j] <- V22[i, j] <- (gamma + delta[i]) * kron(i, j)
      if(i > 1) {
        V11[i, j] <- V11[i, j] - delta[i - 1] * kron(i - 1, j)
        V22[i, j] <- V22[i, j] - delta[i - 1] * kron(i - 1, j)
      }
    }
  }
  F12 <- theta * F11
  F22 <- theta * F21
  
  # Construct the full matrix
  Fmat <- rbind(cbind(F11, F12), cbind(F21, F22))
  Vmat <- rbind(cbind(V11, V12), cbind(V21, V22))
  
  # Next-generation matrix
  if(type == "SIR") {
    G <- F11 %*% solve(V11)
  } else {
    G <- Fmat %*% solve(Vmat)
  }

  return(G)
}