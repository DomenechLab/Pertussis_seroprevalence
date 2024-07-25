#######################################################################################################
# Run and check simulation of pertussis serotransmission model
# All rates are PER YEAR; time unit is YEAR
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
library(pomp)
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
library(reshape2)

# Top-level parameters ----------------------------------------------------
dt_sim <- 1 # Time step separating simulated data points 
dt_mod <- 1e-3 # Time step for stochastic model simulator
country_nm <- c("United-Kingdom", "USA")
n_sims <- 10 # No of stochastic simulations 
n_years_sim <- 300 # No of years of simulations (NB: vaccine is introduced at year 150)
n_years_end <- 20 # No of years to consider at the end of the simulation
rho_V_val <- 0.5 # Probability of immune boosting (from vaccinated state)
alpha_V_val <- 0.02 # Waning rate of vaccine-derived immunity (per year)
vac_cov_prim <- 0.95 # Vaccine coverage from primary series (NB: for boosters, coverage is assumed 10% lower) 
ages_to_vac <- c(2, 3) # Indices of ages to vaccinate

contact_data <- "Mistry" # Source for contact data, either "Mistry" (85 1-yr groups, age 0 to 84) or "Prem" (16 5-yr age groups, age 0-4 to 75-79)
stopifnot(contact_data %in% c("Mistry", "Prem"))
country_nm <- ifelse(contact_data == "Mistry", country_nm[1], country_nm[2])
if(!dir.exists(paste0("_saved/", country_nm))) dir.create(paste0("_saved/", country_nm))
nm_file_save <- sprintf("_saved/%s/alphaV_%.2f-rhoV_%.1f-vacCov_%.2f-%ddoses", 
                        country_nm, alpha_V_val, rho_V_val, vac_cov_prim, length(ages_to_vac))

# Set demographic parameters ----------------------------------------------
# As in Mistry et al., stratification in 1-yr age groups, from age 0 to age 84 (85 age groups overall) 
# Stratify age 0 into two subgroups for the primary vaccination course 
Ntot_val <- 1e7 # Total population size
b_rate <- ifelse(contact_data == "Mistry", 1 / 85, 1 / 80)
nA <- ifelse(contact_data == "Mistry", 86, 81)

delta_vec <- 1 / c(6 / 12, 6 / 12, rep(1, nA - 2)) # Aging rates
N_vec <- b_rate / delta_vec * Ntot_val # Age-specific population sizes
stopifnot(all.equal(sum(N_vec), Ntot_val))

# Data frame with age bounds
age_df <- data.frame(age_fac = 1:nA, age_max = cumsum(1 / delta_vec)) %>% 
  mutate(age_min = c(0, age_max[-length(age_max)]), 
         age_mid = (age_min + age_max) / 2) %>% 
  select(age_fac, age_min, age_mid, age_max)

# Age breaks for aggregated age groups
age_breaks <- c(0, age_df$age_min[2], 1, 5, 10, 20, 40, 60, Inf)

# Add factor age groups: 0-3 mo, 4-11 mo, 1-4 yr, 5-9 yr, 10-19 yr, 20-39 yr, 40-59 yr, >=60 yr
age_df$age_cat <- cut(age_df$age_min, 
                      breaks = age_breaks, 
                      right = F, 
                      include.lowest = T)

# Set contact matrix ----------------------------------------------------------
F_mat <- CreateContactMatrix(country_nm = ifelse(contact_data == "Mistry", country_nm[1], country_nm[2]), 
                             Nvec = N_vec, 
                             source_dat = contact_data, 
                             debug = T)

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
parms[c("rho_R", "alpha_R")] <- c(6.6, 1 / 34)
parms[c("alpha_V", "rho_V")] <- c(alpha_V_val, rho_V_val)
for(i in seq_along(ages_to_vac)) {
  parms[paste0("p_V_", ages_to_vac[i])] <- ifelse(i == 1, vac_cov_prim, vac_cov_prim - 0.1)
}

stopifnot(all(parms >=0))

# Check initial conditions
coef(seroMod, names(parms)) <- unname(parms)
x0 <- rinit(seroMod)
stopifnot(all.equal(sum(x0), Ntot_val))
print(coef(seroMod, c("alpha_V", "rho_V", "alpha_R", "rho_R", "p_V_1", "p_V_2", "p_V_3")))

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
    labs(x = "Time (years)", y = "Age", fill = "Proportion/Rate", title = s)
  print(pl)
}

# Check population sizes --------------------------------------------------
id_check <- 1

pop_checks <- sims_all %>% 
  filter(.id == id_check) %>% 
  select(time, age_mid, pop) %>% 
  unique()

pl <- ggplot(data = pop_checks, 
             mapping = aes(x = time, y = pop / Ntot_val, color = age_mid, group = age_mid)) + 
  geom_line() + 
  theme_classic() + 
  scale_color_viridis(option = "magma", direction = -1) + 
  labs(x = "Time (years)", y = "Count")
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
  labs(x = "Time (years)", y = "Count")
print(pl)

# Calculate mean age at first infection -----------------------------------
MAI_pre_vac <- sims_all %>% 
  filter(var_nm == "Ci1", between(time, parms["t_V"] - 9, parms["t_V"])) %>% 
  group_by(.id, time) %>% 
  summarise(MAI = sum(age_mid * n) / sum(n)) %>% 
  ungroup() %>% 
  group_by(.id) %>% 
  summarise(MAI = mean(MAI)) %>% 
  ungroup()

pl <- ggplot(data = MAI_pre_vac, mapping = aes(x = .id, y = MAI)) + 
  geom_col() + 
  labs(x = "Simulation ID", y = "MAI(prevac) (years)", title = sprintf("Mean(MAI) = %.1f years", mean(MAI_pre_vac$MAI)))
print(pl)

# Plot serological endpoints -------------------------------------------
tmp <- sims_agg %>% 
  filter(time >= max(time) - n_years_end - 1, var_nm %in% c("seroPrev", "seroPPV", "seroInc", "trueInc")) %>% 
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
  labs(x = "Seroprevalence (%)", y = "Age group", title = "Seroprevalence") +
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
  labs(x = "PPV of serology", y = "Age group", title = "PPV") +
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
       title = "Incidence") +
  scale_x_log10() + 
  theme_classic()
print(pl)

#######################################################################################################
# End
#######################################################################################################