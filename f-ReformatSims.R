#######################################################################################################
# Function to reformat simulations from wide format to long format
#######################################################################################################

ReformatSims <- function(df_sim, dt_sim, age_breaks = c(0, 0.1, 1, 5, 10, 20, 40, 60, Inf)) {
  #Args: 
  # df_sim: data frame with simulations in wide format [data frame]
  # dt_sim: time separating observations [numeric]
  # age_breaks: age breaks for broad age groups [numeric vector] 
  # NB: the second value (0.1) is arbitrary, but must be below the age of the second age group   
  # Returns: 
  # Two data frames in long format, one with all age groups, the other with merged age groups
  
  # Compute serological endpoints
  for(i in seq_len(nA)) {
    df_sim[[paste0("seroPrev_", i)]] <- rowSums(df_sim[, paste0(c("Rp1_", "Rp2_", "Vp_"), i)]) # Seroprevalence
    df_sim[[paste0("seroInc_", i)]] <- df_sim[[paste0("Cs_", i)]] / dt_sim # Sero-incidence
    df_sim[[paste0("trueInc_", i)]] <- rowSums(df_sim[, paste0(c("Ci1_", "Ci2_"), i)]) / dt_sim # True incidence
    df_sim[[paste0("seroPPV_", i)]] <- df_sim[[paste0("Rp1_", i)]] / df_sim[[paste0("seroPrev_", i)]] # PPV of serology
  }
  
  # Convert to long format (86 age groups)
  sim_long_all <- df_sim %>% 
    select(-c("X")) %>% 
    pivot_longer(cols = -c(".id", "time"), names_to = "var", values_to = "n") %>% 
    separate(col = "var", into = c("var_nm", "age_fac"), sep = "_", remove = T) %>%
    mutate(var_type = if_else(var_nm %in% state_vars_nm, "state", 
                              if_else(var_nm %in% accum_vars_nm, "accum", "sero")), 
           age_fac = as.integer(age_fac)) %>% 
    left_join(y = age_df) %>% 
    group_by(.id, time, age_fac) %>% 
    mutate(pop = sum(n[var_type == "state"])) %>% 
    ungroup() %>% 
    mutate(age_cat = cut(x = age_min, breaks = age_breaks, include.lowest = T, right = F)) %>% 
    select(.id, time, var_nm, var_type, starts_with("age"), n, everything())
  
  # Aggregate age groups
  sim_long_merge <- sim_long_all %>% 
    group_by(.id, time, var_nm, var_type, age_cat) %>% 
    summarise(n = sum(n), 
              pop = sum(pop)) %>% 
    ungroup()
  
  # Recast in wide format to recalculate seroPPV, then recast in long format
  sim_long_merge <- sim_long_merge %>% 
    select(-var_type) %>% 
    pivot_wider(names_from = c("var_nm"), values_from = "n") %>% 
    mutate(seroPPV = Rp1 / seroPrev) %>% 
    pivot_longer(cols = -c(".id", "time", "age_cat", "pop"), names_to = "var_nm", values_to = "n") %>% 
    mutate(var_type = if_else(var_nm %in% state_vars_nm, "state", 
                              if_else(var_nm %in% accum_vars_nm, "accum", "sero"))) %>% 
    select(.id, time, var_nm, var_type, age_cat, n, pop) %>% 
    arrange(.id, time, var_nm, age_cat)
  
  # Return
  out <- list("all_ages" = sim_long_all, "merged_ages" = sim_long_merge)
  return(out)
  
}