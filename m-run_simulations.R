#######################################################################################################
# Run and check simulation of pertussis serotransmission model
# All rates are PER YEAR; time unit is YEAR
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
library(pomp)
library(reshape2)
source("f-Project_M.R")
source("f-Upsize_M.R")
source("f-PlotMatrix.R")
source("f-CreateSerotransModel.R")
source("f-ReformatSims.R")
source("f-CreateContactMatrix.R")
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
print(packageVersion("pomp"))

# Set fixed parameters ----------------------------------------------------
run_clust <- T

dt_sim <- 1 # Time step separating simulated data points 
dt_mod <- 1e-3 # Time step for stochastic model simulator
n_sims <- 10 # No of stochastic simulations 
n_years_sim <- 300 # No of years of simulations (NB: vaccine is introduced at year 150)
n_years_end <- 20 # No of years to consider at the end of the simulation
ages_to_vac <- c(2, 3, 7) # Indices of ages to vaccinate (age = index - 2)
vac_cov <- c(0.90, 0.90, 0.80) # Age-specific effective vaccine coverage 
stopifnot(length(ages_to_vac) == length(vac_cov))

# SCM data
contact_data <- "Mistry" # Source for contact data, either "Mistry" (85 1-yr groups, age 0 to 84) or "Prem" (16 5-yr age groups, age 0-4 to 75-79)
stopifnot(contact_data %in% c("Mistry", "Prem"))

# Demography
demog_type <- "empirical" # Type of demographic structure: "empirical" (based on actual data) or "type1" (type-I mortality, synthetic population) 
stopifnot(demog_type %in% c("empirical", "type1"))
nA <- 81 # No of age groups
delta_vec_type1 <- 1 / c(2 / 12, 10 / 12, rep(1, nA - 2)) # Aging rates for type I mortality

# Set control parameters --------------------------------------------------
if(run_clust) {
  country_nm <- as.character(Sys.getenv("COUNTRY")); print(country_nm)
  rho_V_val <- as.numeric(Sys.getenv("RHO_V")); print(rho_V_val)
  alpha_V_val <- 1 / as.numeric(Sys.getenv("D_V")); print(alpha_V_val)
  rho_R_val <- as.numeric(Sys.getenv("RHO_R")); print(rho_R_val)
  alpha_R_val <- 1 / as.numeric(Sys.getenv("D_R")); print(alpha_R_val)
  #vac_cov_prim <- as.numeric(Sys.getenv("V_COV")); print(vac_cov_prim)
  vac_cov_prim <- vac_cov[1]
  
} else {
  country_nm <- "Portugal"
  rho_V_val <- 0.25 # Immune boosting coefficient (from V state)
  alpha_V_val <- 0.02 # Waning rate of vaccine-derived immunity (per year)
  rho_R_val <- 6.6 # Immune boosting coefficient (from R state)
  alpha_R_val <- 1 / 34 # Waning rate of infection-derived immunity (per year)
  vac_cov_prim <- vac_cov[1] # Vaccine coverage from primary series (NB: for boosters, coverage is assumed 10% lower) 
}

# Name of file to save
nm_file_save <- sprintf("_outputs_cluster/saved/%s-DV_%.0f-rhoV_%.2f-DR_%.0f-rhoR_%.2f-vacCov_%.2f-%ddoses", 
                        country_nm, 1 / alpha_V_val, rho_V_val, 1 / alpha_R_val, rho_R_val, vac_cov_prim, length(ages_to_vac))
pdf(file = paste0(nm_file_save, ".pdf"), width = 12, height = 8)

# Set demographic parameters, based on real demographic data ---------------------------------------------------
# Stratify age 0 into two subgroups for the primary vaccination course 

if(demog_type == "empirical") {
  demog_dat <- read_csv(file = sprintf("_data/_demog/_2010/%s_country_level_age_distribution_85.csv", country_nm), 
                        col_names = c("age", "pop"), 
                        col_types = "d") %>% 
    arrange(age) %>% 
    filter(age <= 79)
  
  # Population sizes
  Ntot_val <- sum(demog_dat$pop) # Total population size
  N_vec <- c(demog_dat$pop[1] / delta_vec_type1[1], demog_dat$pop[1] / delta_vec_type1[2],  demog_dat$pop[-1]) # Age-specific population sizes
  
  # Birth rate
  b_rate <- read_csv2(file = "_data/_demog/birth_rates_2010.csv", col_names = T, col_types = "cd")
  b_rate <- b_rate$birth_rate[b_rate$country == country_nm] / 1e3
  stopifnot(length(b_rate) == 1 && is.numeric(b_rate))
  
  # Aging rates
  delta_vec <- numeric(nA)
  delta_vec[1] <- b_rate * Ntot_val / N_vec[1]
  delta_vec[-1] <- delta_vec[1] * N_vec[1] / N_vec[-1]
  
  stopifnot(length(N_vec) == nA)
  stopifnot(all.equal(sum(N_vec), Ntot_val))
}

# Set demographic parameters, type I mortality in synthetic populations ----------------------------------------------
if(demog_type == "type1") {
  Ntot_val <- 1e7 # Total population size
  b_rate <- 1 / (nA - 1) # Birth rate (per year)
  
  delta_vec <- delta_vec_type1 # Aging rates
  N_vec <- b_rate / delta_vec * Ntot_val # Age-specific population sizes
  stopifnot(all.equal(sum(N_vec), Ntot_val))
  demog_dat <- data.frame(age = 0:(nA - 2), pop = c(sum(N_vec[1:2]), N_vec[-c(1:2)]))
}

# Define age groups ----------------------------------------------
# Data frame with age bounds
age_df <- data.frame(age_fac = 1:nA, age_max = cumsum(1 / delta_vec_type1)) %>% 
  mutate(age_min = c(0, age_max[-length(age_max)]), 
         age_mid = (age_min + age_max) / 2) %>% 
  select(age_fac, age_min, age_mid, age_max)

# Age breaks for aggregated age groups
age_breaks <- c(0, age_df$age_min[2], 1, 5, 10, 15, 20, 25, 40, 50, 60, Inf)

# Add factor age groups: 0-3 mo, 4-11 mo, 1-4 yr, 5-9 yr, 10-19 yr, 20-39 yr, 40-59 yr, >=60 yr
age_df$age_cat <- cut(age_df$age_min, 
                      breaks = age_breaks, 
                      right = F, 
                      include.lowest = T)

# Set contact matrix ----------------------------------------------------------
SCMs <- CreateContactMatrix(country_nm = country_nm, 
                            Nvec = N_vec, 
                            source_dat = contact_data, 
                            trim_mat = T,  
                            debug = F)
F_mat <- SCMs$F_mat

# Create and initialize POMP model -------------------------------------------------------
seroMod <- CreateSerotransMod(nA = nA, 
                              f_mat = F_mat, 
                              dt_mod = dt_mod, 
                              dt_sim = dt_sim, 
                              t_sim = n_years_sim, 
                              debug_bool = F)
parms <- coef(seroMod)

# Set initial conditions (assume full susceptibility at time 0s)
N_vec <- round(N_vec)
Ntot_val <- sum(N_vec)
e0_vec <- rep(100L, nA)
s0_vec <- as.integer(N_vec - e0_vec)
parms[str_detect(names(parms), "S1_[0-9]+_0")] <- s0_vec
parms[str_detect(names(parms), "E1_[0-9]+_0")] <- e0_vec

# Set demographic parameters
parms[paste0("delta_", 1:nA)] <- delta_vec
parms[paste0("N_", 1:nA)] <- N_vec
parms["N_tot"] <- Ntot_val
parms["b_rate"] <- b_rate
parms[c("alpha_R", "rho_R")] <- c(alpha_R_val, rho_R_val)
parms[c("alpha_V", "rho_V")] <- c(alpha_V_val, rho_V_val)
parms[paste0("p_V_", ages_to_vac)] <- vac_cov

stopifnot(all(parms >=0))

# Check initial conditions
coef(seroMod, names(parms)) <- unname(parms)
x0 <- rinit(seroMod)
stopifnot(all.equal(sum(x0), Ntot_val))
print(sprintf("Country: %s, Ntot = %.1fM, b_rate=%.3f per yr", country_nm, Ntot_val / 1e6, b_rate))
print(coef(seroMod, c("alpha_V", "rho_V", "alpha_R", "rho_R", paste0("p_V_", ages_to_vac))))

# Names of state variables
state_vars_nm <- c("S1", "S2", "E1", "E2", "I1", 
                   "I2", "R", "Re", "Rp1", "Rp2", 
                   "V", "Ve", "Vp")
accum_vars_nm <- c("Ci1", "Ci2", "Cs")

# Run simulation -----------------------------------------------------
sims_list <- bake(file = paste0(nm_file_save, ".rds"), 
                  seed = 2186L, 
                  expr = {
                    simulate(seroMod, nsim = n_sims, format = "data.frame") %>% 
                      ReformatSims(dt_sim = dt_sim, age_breaks = age_breaks)
                  })

vars_nm <- unique(sims_list[[1]]$var_nm)
sims_all <- sims_list[[1]]
sims_agg <- sims_list[[2]]

# Check convergence -------------------------------------------------------
vars_checks <- c("S1", "S2", "trueInc", "seroInc", "seroPrev")
id_check <- 1

for(s in vars_checks) {
  pl <- ggplot(data = sims_all %>% filter(var_nm == s, .id == id_check, time >= 50), 
               mapping = aes(x = time, y = age_fac, fill = n / pop)) + 
    geom_tile() + 
    scale_fill_viridis(option = "magma", direction = -1) + 
    labs(x = "Time (years)", y = "Age", fill = "Proportion/Rate", title = paste0(country_nm, ", ", s))
  print(pl)
}

# Check population sizes --------------------------------------------------
id_check <- 1

pop_checks <- sims_all %>% 
  filter(.id == id_check) %>% 
  select(time, age_min, pop) %>% 
  unique()

# Population over time
pl <- ggplot(data = pop_checks, 
             mapping = aes(x = time, y = pop / Ntot_val, color = age_min, group = age_min)) + 
  geom_line() + 
  theme_classic() + 
  scale_color_viridis(option = "magma", direction = -1) + 
  labs(x = "Time (years)", y = "Count", title = paste0(country_nm, ", Population over time (all age groups)"))
print(pl)

# Population pyramid
pl <- ggplot(data = pop_checks, mapping = aes(x = age_min, y = pop / 1e6, group = age_min)) + 
  geom_boxplot() + 
  geom_point(data = demog_dat, 
             mapping = aes(x = age, y = pop / 1e6, group = NULL), color = "red") + 
  labs(x = "Age (years)", y = "Population size (millions)", title = paste0(country_nm, ", Population pyramid (red: data)"))
print(pl)

pop_checks <- sims_agg %>% 
  filter(.id == id_check) %>% 
  select(time, age_cat, pop) %>% 
  unique()

pl <- ggplot(data = pop_checks, 
             mapping = aes(x = time, y = pop / Ntot_val, color = age_cat, group = age_cat)) + 
  geom_line() + 
  theme_classic() + 
  scale_color_viridis(option = "magma", direction = -1, discrete = T) + 
  labs(x = "Time (years)", y = "Count", title = paste0(country_nm, ", Population over time (aggregated age groups)"))
print(pl)

# Calculate mean age at first infection (pre-vaccine era) -----------------------------------

# Calculate age distribution of primary and primary infections
# Average over time: inc = E_t(cases_t / pop_t)
age_dist <- sims_all %>% 
  filter(var_nm %in% c("Ci1", "Ci2"), 
         between(time, parms["t_V"] - n_years_end + 1, parms["t_V"])) %>% 
  group_by(.id, var_nm, age_min, age_mid, age_max) %>% 
  summarise(inc_mean = mean(n / pop), 
            inc_sd = sd(n / pop)) %>% 
  ungroup()

# Calculate mean and mode of age distribution
age_dist_stats <- age_dist %>% 
  filter(var_nm == "Ci1") %>% 
  group_by(.id) %>% 
  summarise(mean = sum(age_mid * inc_mean) / sum(inc_mean), 
            mode = age_mid[which.max(inc_mean)]) %>% 
  ungroup()

# Plot age distribution
# Add mean (SD) for mean and mode across simulations
pl <- ggplot(data = age_dist, 
             mapping = aes(x = age_min, y = 1e2 * inc_mean, group = interaction(.id, var_nm), color = var_nm)) + 
  geom_line() + 
  scale_y_sqrt() + 
  theme_classic() + 
  labs(x = "Age (years)", 
       y = "Incidence rate (per year per 100)", 
       title = sprintf("%s, Prevaccine era, mean(A1) = %.1f (%.1f), mode(A1) = %.1f (%.1f) years", 
                       country_nm, 
                       mean(age_dist_stats$mean), sd(age_dist_stats$mean), 
                       mean(age_dist_stats$mode), sd(age_dist_stats$mode)))
print(pl)

# Calculate mean age at second infection (vaccine era) -----------------------------------

# Calculate age distribution of primary and primary infections
# Average over time: inc = E_t(cases_t / pop_t)
age_dist <- sims_all %>% 
  filter(var_nm %in% c("Ci1", "Ci2"), 
         between(time, max(time) - n_years_end + 1, max(time))) %>% 
  group_by(.id, var_nm, age_min, age_mid, age_max) %>% 
  summarise(inc_mean = mean(n / pop), 
            inc_sd = sd(n / pop)) %>% 
  ungroup()

# Calculate mean and mode of age distribution
age_dist_stats <- age_dist %>% 
  filter(var_nm == "Ci2") %>% 
  group_by(.id) %>% 
  summarise(mean = sum(age_mid * inc_mean) / sum(inc_mean), 
            mode = age_mid[which.max(inc_mean)]) %>% 
  ungroup()

# Plot age distribution
# Add mean (SD) for mean and mode across simulations
pl <- ggplot(data = age_dist, 
             mapping = aes(x = age_min, y = 1e5 * inc_mean, group = interaction(.id, var_nm), color = var_nm)) + 
  geom_line() + 
  scale_y_sqrt() + 
  theme_classic() + 
  labs(x = "Age (years)", 
       y = "Incidence rate (per year per 100,000)", 
       title = sprintf("%s, Vaccine era, mean(A2) = %.1f (%.1f), mode(A2) = %.1f (%.1f) years", 
                       country_nm, 
                       mean(age_dist_stats$mean), sd(age_dist_stats$mean), 
                       mean(age_dist_stats$mode), sd(age_dist_stats$mode)))
print(pl)

# Plot serological endpoints -------------------------------------------
tmp <- sims_agg %>% 
  filter(between(time, max(time) - n_years_end + 1, max(time)), 
         var_nm %in% c("seroPrev", "seroPPV", "seroInc", "trueInc")) %>% 
  select(-var_type) %>% 
  mutate(prop = n / pop)

# Summary: median
tmp_sumry <- tmp %>% 
  group_by(var_nm, age_cat) %>% 
  summarise(med_n = median(n), 
            med_prop = median(prop)) %>% 
  ungroup()

# Seroprevalence
pl <- ggplot(data = tmp %>% filter(var_nm == "seroPrev"), 
             mapping = aes(x = 1e2 * prop, y = age_cat)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      #stat = "density",
                      jittered_points = T,
                      alpha = 0.5,
                      point_shape = "|", 
                      position = position_points_jitter(height = 0)) + 
  geom_text(mapping = aes(x = 1e2 * med_prop, y = age_cat, label = round(1e2 * med_prop, 1)), 
            data = tmp_sumry %>% filter(var_nm == "seroPrev"), 
            color = "red") + 
  labs(x = "Seroprevalence (%)", y = "Age group", title = paste0(country_nm, ", Seroprevalence")) +
  theme_classic()
print(pl)

# PPV of serology
pl <- ggplot(data = tmp %>% filter(var_nm == "seroPPV"), 
             mapping = aes(x = 1e2 * n, y = age_cat)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      alpha = 0.5,
                      jittered_points = T, 
                      point_shape = "|", 
                      position = position_points_jitter(height = 0)) + 
  geom_text(mapping = aes(x = 1e2 * med_n, y = age_cat, label = round(1e2 * med_n, 1)), 
            data = tmp_sumry %>% filter(var_nm == "seroPPV"), 
            color = "red") + 
  labs(x = "PPV of serology", y = "Age group", title = paste0(country_nm, ", PPV")) +
  theme_classic()
print(pl)

# Sero-incidence vs. true incidence
pl <- ggplot(data = tmp %>% filter(var_nm %in% c("seroInc", "trueInc")), 
             mapping = aes(x = 1e5 * prop, y = age_cat, color = var_nm, fill = var_nm)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      jittered_points = F,
                      stat = "density_ridges", 
                      scale = 0.9,
                      alpha = 0.5) + 
  geom_text(mapping = aes(x = 1e5 * med_prop, y = age_cat, label = round(1e5 * med_prop), color = var_nm), 
            data = tmp_sumry %>% filter(var_nm %in% c("seroInc", "trueInc"))) + 
  labs(x = "Incidence rate (per year per 100,000)", 
       y = "Age group",
       title = paste0(country_nm, ", Incidence")) +
  scale_x_log10() + 
  theme_classic()
print(pl)

# Exit instructions -------------------------------------------------------
dev.off()

#######################################################################################################
# End
#######################################################################################################