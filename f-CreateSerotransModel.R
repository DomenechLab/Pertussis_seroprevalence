#######################################################################################################
# Create pomp object for the pertussis serotransmission model
# The time unit is YEARS; all rates are PER YEAR
#######################################################################################################

CreateSerotransMod <- function(nA, 
                               f_mat,
                               dt_mod = 0.01, 
                               dt_sim = 0.1,
                               t_sim = 200, 
                               debug_bool = F) {
  # Args: 
  # nA: no of age classes [integer]
  # f_mat: contact matrix (density scale) [matrix]
  # dt_mod: time step (in years) for stochastic model simulator [numeric]
  # dt_sim: time step (in years) for simulated data [numeric]
  # t_sim: no of years of simulation [numeric]
  # debug_bool: should messages be displayed to help debug? [boolean] 
  
  # Check contact matrix
  stopifnot(nrow(f_mat) == ncol(f_mat))
  stopifnot(nrow(f_mat) == nA)
  stopifnot(isSymmetric(f_mat))
  
  # Extract C code from file
  mod_code <- readLines("c-serotransmission_model.c")
  components_nm <- c("globs", "stochSim")
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for(nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0("start_", nm)):str_which(mod_code, paste0("end_", nm))] %>% 
      str_flatten(collapse = "\n")
    
    if(nm == "globs") {
      components_l[[nm]] <- paste(components_l[[nm]], 
                                  sprintf("static int nA = %d;\nstatic int debug = %d;", 
                                          nA, as.integer(debug_bool)), 
                                  sep = "\n")
    }
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  # Names of model variables 
  states_nm <- c("S1", "S2", "E1", "E2", "I1", 
                 "I2", "R", "Re", "Rp1", "Rp2", 
                 "V", "Ve", "Vp", "Ci1", "Ci2", 
                 "Cs")
  all_states_nm <- sapply(paste0(states_nm, "_"), paste0, 1:nA) %>% as.vector()
  
  # Names of parameters
  # NB: negative values indicate parameters that have no default value and must be assigned
  params_vec <- c(
    setNames(object = as.numeric(t(f_mat)), nm = paste0("f_", 1:(nA * nA))), # Contact matrix
    setNames(object = rep(0, nA), paste0("p_V_", 1:nA)), # Effective vaccination coverage (proportion)
    setNames(object = rep(-1, nA), paste0("delta_", 1:nA)), # Aging rates
    setNames(object = rep(-1, nA), paste0("N_", 1:nA)), # Age-specific populations
    "N_tot" = -1, # Total population size
    "b_rate" = -1, # Birth rate (per capita)
    "t_V" = 150, # Start time of vaccination
    "alpha_V" = 0.01067172, # Rate of waning vaccine-derived immunity
    "q1" = 0.09373881, # Susceptibility in 0-9 yo
    "q21" = 0.5117505, # Susceptibility in 10-19 yo relative to that in 0-9 yo
    "q32" = 0.1725993 , # Susceptibility in >=20 yo relative to that in 10-19 yo
    "theta" = 0.9874919, # Transmissibility of secondary infections relative to that of primary infections
    "iota" = 1e-3, # Prevalence of imported cases
    "sigma" = 365 / 8, # Latent rate 
    "gamma" = 365 / 15, # Recovery rate
    "alpha_R" = 1 / 34, # Rate of waning infection-derived immunity (Table 1, Lavine et al., 2013)
    "rho_R" = 6.6, # Boosting coefficient (Table 1, Lavine et al., 2013)
    "rho_V" = 1, # Boosting coefficient
    "t_p_R" = 23 / 365, # Average time to seroconversion after immune boost
    "t_p_V" = 23 / 365, # # Average time to seroconversion after immune boost
    "t_n_I" = 1, # Average duration of seropositivity after infection
    "t_n_R" = 1, # Average duration of seropositivity after immune boost from R state
    "t_n_V" = 1, # Average duration of seropositivity after immune boost from V state
    setNames(object = rep(0, length(all_states_nm)), nm = paste0(all_states_nm, "_0")) # Initial conditions
  )
  
  #browser()
  po <- pomp(data = data.frame(time = seq(from = dt_sim, to = t_sim, by = dt_sim), 
                               X = NA),
             times = "time",
             t0 = 0,
             obsnames = "X",
             statenames = all_states_nm,
             accumvars = c(
               paste0("Ci1_", 1:nA),
               paste0("Ci2_", 1:nA), 
               paste0("Cs_", 1:nA)
             ),
             paramnames = names(params_vec), 
             params = params_vec,
             globals = components_l[["globs"]],
             rprocess = euler(step.fun = components_l[["stochSim"]], delta.t = dt_mod), 
             rinit = NULL
             #cdir = "_help/", 
             #cfile = "codes"
  )
  return(po)
}

#######################################################################################################
# END
#######################################################################################################
