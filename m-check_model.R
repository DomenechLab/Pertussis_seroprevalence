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
source("f-compute_R0.R")
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
print(packageVersion("pomp"))
library(reshape2)

# Top-level parameters ----------------------------------------------------
dt_sim <- 1 # Time step separating simulated data points 
dt_mod <- 1e-3 # Time step for stochastic model simulator
#country_nm <- c("Israel", "USA")
country_nm <- Sys.getenv("country_nm")
n_sims <- 1 # No of stochastic simulations 
rho_V_val <- 0.25 # Probability of immune boosting
alpha_V_val <- 0.02 # Waning rate of infection-derived immunity
contact_data <- "Mistry" # Source for contact data, either "Mistry" (85 1-yr groups, age 0 to 84) or "Prem" (16 5-yr age groups, age 0-4 to 75-79)
stopifnot(contact_data %in% c("Mistry", "Prem"))
print(country_nm)

# Set demographic parameters ----------------------------------------------
# As in Mistry et al., stratification in 1-yr age groups, from age 0 to age 84 (85 age groups overall) 
# Stratify age 0 into two subgroups for the primary vaccination course 
Ntot_val <- 1e7 # Total population size
b_rate <- ifelse(contact_data == "Mistry", 1 / 85, 1 / 80)
nA <- ifelse(contact_data == "Mistry", 86, 81)

delta_vec <- 1 / c(6 / 12, 6 / 12, rep(1, nA - 2)) # Aging rates
N_vec <- b_rate / delta_vec * Ntot_val # Age-specific population sizes
stopifnot(all.equal(sum(N_vec), Ntot_val))

# Data frame wti age bounds
age_df <- data.frame(age_fac = 1:nA, age_max = cumsum(1 / delta_vec)) %>% 
  mutate(age_min = c(0, age_max[-length(age_max)]), 
         age_mid = (age_min + age_max) / 2) %>% 
  select(age_fac, age_min, age_mid, age_max)

age_breaks <- c(0, age_df$age_min[2], 1, 5, 10, 20, 40, 60, Inf)

# Add factor age groups: 0-3 mo, 4-11 mo, 1-4 yr, 5-9 yr, 10-19 yr, 20-39 yr, 40-59 yr, >=60 yr
age_df$age_cat <- cut(age_df$age_min, 
                      breaks = age_breaks, 
                      right = F, 
                      include.lowest = T)

# Set contact matrix ----------------------------------------------------------
SCMs <- CreateContactMatrix(country_nm = ifelse(contact_data == "Mistry", country_nm[1], country_nm[2]), 
                             Nvec = N_vec, 
                             source_dat = contact_data, 
                             debug = T)
F_mat <- SCMs$F_mat

# Create and initialize POMP model -------------------------------------------------------
seroMod <- CreateSerotransMod(nA = nA, 
                              f_mat = F_mat, 
                              dt_mod = dt_mod, 
                              dt_sim = dt_sim, 
                              t_sim = 300, 
                              debug_bool = F)
parms <- coef(seroMod)

# Set initial conditions
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
stopifnot(all(parms >=0))

# Check initial conditions
coef(seroMod, names(parms)) <- unname(parms)
x0 <- rinit(seroMod)
stopifnot(all.equal(sum(x0), Ntot_val))

# Names of state variables
state_vars_nm <- c("S1", "S2", "E1", "E2", "I1", 
                   "I2", "R", "Re", "Rp1", "Rp2", 
                   "V", "Ve", "Vp")
accum_vars_nm <- c("Ci1", "Ci2", "Cs")

# Run test simulation -----------------------------------------------------
vac_cov <- 0.95
coef(seroMod, names(parms)) <- unname(parms)
coef(seroMod, c("rho_R", "alpha_R")) <- c(6.6, 1 / 34)
coef(seroMod, c("p_V_2", "p_V_3")) <- c(vac_cov, vac_cov - 0.1) # Vaccinate second and third age group
coef(seroMod, c("alpha_V", "rho_V")) <- c(alpha_V_val, rho_V_val) # Vaccinate second age group

sim_test <- simulate(seroMod, nsim = n_sims, format = "data.frame")

sims_list <- ReformatSims(df_sim = sim_test, dt_sim = dt_sim, age_breaks = age_breaks)
vars_nm <- unique(sims_list[[1]]$var_nm)

# Calculate mean age at first infection and R0 -----------------------------------
id_cur <- 1
sim_cur <- filter(sims_list[[1]], .id == id_cur)

q1 <- unname(coef(seroMod, "q1"))
q2 <- unname(coef(seroMod, "q21") * q1)
q3 <- unname(coef(seroMod, "q32") * q2)

q_vec <- c(rep(q1, 11), rep(q2, 10), rep(q3, nA - 21))

NGM <- compute_R0(theta = unname(coef(seroMod, "theta")), 
                  gamma = unname(coef(seroMod, "gamma")), 
                  N = N_vec, 
                  q = q_vec, 
                  delta = delta_vec, 
                  Cmat = SCMs$M_mat)

R0_val <- NGM %>% eigen() %>% pluck("values") %>% abs() %>% max()

MAI <- sim_cur %>% 
  filter(var_nm == "Ci1") %>% 
  group_by(time) %>% 
  summarise(MAI = sum(age_mid * n) / sum(n)) %>% 
  ungroup() %>% 
  filter(time >= 50)

MAI_pre_vac <- MAI %>% 
  filter(between(time, parms["t_V"] - 9, parms["t_V"])) %>% 
  pluck("MAI") %>% 
  mean()

plot(MAI$time, MAI$MAI, type = "l", 
     xlab = "Time (years)", 
     ylab = "Mean age of first infection (years)", 
     main = sprintf("Country = %s, MAI(prevaccine era) = %.2f years, R0=%.1f", country_nm, MAI_pre_vac, R0_val))
abline(h = MAI_pre_vac, col = "grey")

# Check population sizes for aggregated age groups ------------------------
df_cur <- sims_list[["merged_ages"]] %>% 
  filter(.id == id_cur) %>%  
  select(time, age_cat, pop) %>% 
  unique()

# Population size (age-specific)
pl <- ggplot(data = df_cur, 
             mapping = aes(x = time, y = pop, color = age_cat, group = age_cat)) + 
  geom_line() + 
  theme_classic() + 
  scale_color_viridis(option = "magma", direction = -1, discrete = T) + 
  labs(x = "Time (years)", y = "Count") + 
  ggtitle("Population sizes (age-specific)")
print(pl)

df_cur2 <- df_cur %>% 
  group_by(time) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()

# Total population size
pl <- ggplot(data = df_cur2, 
             mapping = aes(x = time, y = pop / Ntot_val)) + 
  geom_line() + 
  theme_classic() + 
  labs(x = "Time (years)", y = "Count") + 
  ggtitle("Population size (total)")
print(pl)

# Check population sizes for all age groups -------------------------------------------------------------
df_cur <- sims_list[["all_ages"]] %>% 
  filter(.id == id_cur) %>%  
  select(time, age_mid, pop) %>% 
  unique()

# Population size (age-specific)
pl <- ggplot(data = df_cur, 
             mapping = aes(x = time, y = pop, color = age_mid, group = age_mid)) + 
  geom_line() + 
  theme_classic() + 
  scale_color_viridis(option = "magma", direction = -1) + 
  labs(x = "Time (years)", y = "Count") + 
  ggtitle("Population sizes (age-specific)")
print(pl)

df_cur2 <- df_cur %>% 
  group_by(time) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()

# Total population size
pl <- ggplot(data = df_cur2, 
             mapping = aes(x = time, y = pop / Ntot_val)) + 
  geom_line() + 
  theme_classic() + 
  labs(x = "Time (years)", y = "Count") + 
  ggtitle("Population size (total)")
print(pl)

# Susceptibility profile --------------------------------------------------
sim_cur <- sims_list[["all_ages"]] %>% 
  filter(.id == id_cur)

for(s in c("S1", "S2")) {
  pl <- ggplot(data = sim_cur %>% filter(var_nm == s), 
               mapping = aes(x = time, y = age_fac, fill = n / pop)) + 
    geom_tile() + 
    scale_fill_viridis(option = "magma", direction = -1) + 
    labs(x = "Time (years)", y = "Age", fill = "Proportion", title = s)
  print(pl)
}

# Plot all variables (all age groups)------------------------------------------------------
sim_cur <- sims_list[["all_ages"]] %>% 
  filter(.id == id_cur)
t_range <- c(281, 300) 
var_type <- unique(sim_cur$var_type)

for(v in var_type) {
  df_cur <- sim_cur %>% filter(var_type == v)
  var_nm_cur <- unique(df_cur$var_nm)
  
  for(s in var_nm_cur) {
    
    if(s == "seroPPV") {
      pl <- ggplot(data = sim_cur %>% 
                     filter(var_nm == s, between(time, t_range[1], t_range[2])), 
                   mapping = aes(x = time, y = n, color = age_mid, group = age_mid))
    } else {
      pl <- ggplot(data = sim_cur %>% 
                     filter(var_nm == s, between(time, t_range[1], t_range[2])), 
                   mapping = aes(x = time, y = n / pop, color = age_mid, group = age_mid))
    }
    pl <- pl + 
      geom_line() + 
      theme_classic() + 
      scale_color_viridis(option = "magma", direction = -1) + 
      scale_y_sqrt() + 
      labs(x = "Time (years)", y = "Proportion", color = "Age") + 
      ggtitle(sprintf("Variable type: %s, variable name: %s", v, s))
    if(s == "V") {
      pl <- pl + geom_hline(yintercept = vac_cov, color = "grey")
    }
    print(pl)
  }
}

# Plot serological endpoints --------------------------------------------
t_range <- c(281, 300) 
tmp <- sims_list[["all_ages"]] %>% 
  filter(.id == id_cur) %>% 
  filter(var_nm %in% c("seroPrev", "seroPPV"), between(time, t_range[1], t_range[2]))

pl <- ggplot(data = tmp %>% filter(var_nm == "seroPrev"), 
             mapping = aes(x = age_mid, y = 1e2 * n / pop, group = age_mid)) + 
  geom_boxplot(outlier.colour = "grey") + 
  theme_classic() + 
  scale_x_sqrt() + 
  labs(x = "Age (years)", y = "Seroprevalence (%)", title = "Seroprevalence (all age groups)")
print(pl)

pl <- ggplot(data = tmp %>% filter(var_nm == "seroPPV"), 
             mapping = aes(x = age_mid, y = 1e2 * n, group = age_mid)) + 
  geom_boxplot(outlier.colour = "grey") + 
  theme_classic() + 
  scale_x_sqrt() + 
  labs(x = "Age (years)", y = "PPV of serology (%)", title = "PPV (all age groups)")
print(pl)

# Plot all variables (merged age groups)------------------------------------------------------
t_range <- c(281, 300) 

sim_cur <- sims_list[["merged_ages"]] %>% 
  filter(.id == id_cur, between(time, t_range[1], t_range[2]))

var_type <- unique(sim_cur$var_type)

for(v in var_type) {
  df_cur <- sim_cur %>% filter(var_type == v)
  var_nm_cur <- unique(df_cur$var_nm)
  
  for(s in var_nm_cur) {
    
    if(s == "seroPPV") {
      pl <- ggplot(data = sim_cur %>% filter(var_nm == s), 
                   mapping = aes(x = time, y = n, color = age_cat, group = age_cat))
      
      pl2 <- ggplot(data = sim_cur %>% filter(var_nm == s), 
                    mapping = aes(x = n, y = age_cat)) 
      
    } else {
      pl <- ggplot(data = sim_cur %>% 
                     filter(var_nm == s), 
                   mapping = aes(x = time, y = n / pop, color = age_cat, group = age_cat))
      
      pl2 <- ggplot(data = sim_cur %>% filter(var_nm == s), 
                    mapping = aes(x = n / pop, y = age_cat)) 
    }
    pl <- pl + 
      geom_line() + 
      theme_classic() + 
      scale_color_viridis(option = "magma", direction = -1, discrete = T) + 
      scale_y_sqrt() + 
      labs(x = "Time (years)", y = "Proportion", color = "Age") + 
      ggtitle(sprintf("Variable type: %s, variable name: %s", v, s))
    
    pl2 <- pl2 + 
      geom_density_ridges(quantile_lines = T, 
                          quantiles = 2, 
                          jittered_points = T, 
                          alpha = 0.5,
                          position = position_points_jitter(height = 0)) + 
      theme_classic() + 
      labs(x = "Proportion", y = "Age (years)") + 
      ggtitle(sprintf("Variable type: %s, variable name: %s", v, s))
    
    if(s == "V") {
      pl <- pl + geom_hline(yintercept = vac_cov, color = "grey")
    }
    print(pl)
    print(pl2)
  }
}

# Plot contributions to seroprevalence ----------------------------------
tmp <- sims_list[["merged_ages"]] %>% 
  filter(.id == id_cur, 
         between(time, t_range[1], t_range[2]), 
         var_nm %in% c("seroPrev", "Rp1", "Rp2", "Vp", "S1", "S2", "R", "V")) %>% 
  select(-var_type) %>% 
  mutate(prop = n / pop) %>% 
  select(-c(pop, n)) %>% 
  pivot_wider(names_from = "var_nm", values_from = "prop") %>% 
  mutate(Rp1 = Rp1 / seroPrev, 
         Rp2 = Rp2 / seroPrev, 
         Vp = Vp / seroPrev, 
         S = S1 + S2) %>% 
  pivot_longer(cols =-c(".id", "time", "age_cat")) %>% 
  group_by(age_cat, name) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

# Variables S, V, and R
pl <- ggplot(data = tmp %>% filter(name %in% c("S", "R", "V")) %>% mutate(name = fct_relevel(name, "R", after = 1)), 
             mapping = aes(x = age_cat, y = value)) + 
  geom_col(position = "dodge") + 
  facet_wrap(~ name, scales = "fixed", ncol = 2) + 
  labs(x = "Age", y = "Prevalence")
print(pl)

# Variables S, V, and R
pl <- ggplot(data = tmp %>% filter(name %in% c("Rp1", "Rp2", "Vp")), 
             mapping = aes(x = age_cat, y = value)) + 
  geom_col(position = "dodge") + 
  facet_wrap(~ name, scales = "fixed", ncol = 2) + 
  labs(x = "Age", y = "Proportion (relative to seroprevalence)")
print(pl)

# Plot serological endpoints ----------------------------------------------
tmp <- sims_list[["merged_ages"]] %>% 
  filter(.id == id_cur, 
         between(time, t_range[1], t_range[2]), 
         var_nm %in% c("seroPrev", "seroPPV", "seroInc", "trueInc", "Rp1", "Rp2", "Vp")) %>% 
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
                      scale = 0.9, 
                      jittered_points = F, 
                      #point_shape = "|", 
                      position = position_points_jitter(height = 0)) + 
  geom_text(mapping = aes(x = 1e2 * med_prop, y = age_cat, label = round(1e2 * med_prop, 1)), 
            data = tmp_sumry %>% filter(var_nm == "seroPrev")) + 
  labs(x = "Seroprevalence", y = "Age group", title = "Seroprevalence (merged age groups)") +
  theme_classic()
print(pl)

# PPV of serology
pl <- ggplot(data = tmp %>% filter(var_nm == "seroPPV"), 
             mapping = aes(x = 1e2 * n, y = age_cat)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      scale = 0.9,
                      jittered_points = F, 
                      #point_shape = "|", 
                      position = position_points_jitter(height = 0)) + 
  geom_text(mapping = aes(x = 1e2 * med_n, y = age_cat, label = round(1e2 * med_n, 1)), 
            data = tmp_sumry %>% filter(var_nm == "seroPPV")) + 
  labs(x = "PPV of serology", y = "Age group", title = "PPV (merged age groups)") +
  theme_classic()
print(pl)

# Sero-incidence vs. true incidence
pl <- ggplot(data = tmp %>% filter(var_nm %in% c("seroInc", "trueInc")), 
             mapping = aes(x = 1e5 * prop, y = age_cat, color = var_nm, fill = var_nm)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      jittered_points = F,
                      scale = 0.9,
                      alpha = 0.5) + 
  geom_text(mapping = aes(x = 1e5 * med_prop, y = age_cat, label = round(1e5 * med_prop), color = var_nm), 
            data = tmp_sumry %>% filter(var_nm %in% c("seroInc", "trueInc"))) + 
  labs(x = "Incidence (per year per 100,000)", 
       y = "Age group",
       title = "Incidence (merged age groups)") +
  scale_x_log10() + 
  theme_classic()
print(pl)

#######################################################################################################
# End
#######################################################################################################