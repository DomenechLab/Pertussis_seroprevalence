#######################################################################################################
# Run simulations of alternative model with two waning stages and two boosting parameters
# For simplicity, the model is assumed homogeneous (no age structure)
#######################################################################################################

# Load packages -----------------------------------------------------------
rm(list = ls())
source("s-base_packages.R")
library(pomp)
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
print(packageVersion("pomp"))

# Model 1: Base model, homogeneous version --------------------------------------------------------------
# For simplicity, exposed states (Ve and Re) and secondary infections (E2 and I2) are omitted

# Name of state variables
mod1_s_nm <- c("S", "E", "I", "Rp1", "R", "Rp2", "V", "Vp")

# Deterministic skeleton (ODEs)
mod1_skel <- Csnippet("
    // Force of infection
    double lambda = Rzero * (gamma + mu) * I; // Force of infection
    double v_rate = pV / (1.0 - pV) * mu; // Effective vaccine coverage
    
    // ODEs
    DS = mu + alpha * R + alpha * V - (lambda + v_rate + mu) * S;
    DE = lambda * S - (sigma + mu) * E;
    DI = sigma * E - (gamma + mu) * I;
    DRp1 = gamma * I - (1.0 / tn + mu) * Rp1;
    DR = Rp1 / tn + Rp2 / tn - (rho * lambda + alpha + mu) * R;
    DRp2 = rho * lambda * R - (1.0 / tn + mu) * Rp2;
    DV = v_rate * S + Vp / tn - (rho * lambda + alpha + mu) * V;
    DVp = rho * lambda * V - (1.0 / tn + mu) * Vp;
  ")

# Initializer (based on equilibria of VSIRS model with vaccination v_rate)
mod1_init <- Csnippet("
  S = 1.0 / Rzero; 
  V = pV / (1.0 - pV) * mu / (alpha + mu) / Rzero; 
  E = 0.0; 
  I = (1.0 - S - V) / (1.0 + gamma / (alpha + mu));
  R = 1.0 - V - S - I; 
  Rp1 = 0.0; 
  Rp2 = 0.0; 
  Vp = 0.0; 
")

# Base model parameters
parms_mod1 <- c("mu" = 1 / 80, # Birth/death rate
                "pV" = 0.9, # Effective vaccination coverage
                "Rzero" = 15, # Basic reproduction number
                "sigma" = 365 / 8, # 1 / average latent period
                "gamma" = 365 / 15, # 1 / average infectious period
                "tn" = 1.9, # Average duration of seropositivity
                "alpha" = 1 / 40, # 1 / average duration of immunity
                "rho" = 0.5) # Boosting coefficient

# Create model
mod1 <- pomp(
  data = data.frame(time = seq(1, 5e2, by = 1), cases = NA),
  obsnames = "cases", 
  times = "time",
  t0 = 0, 
  rinit = mod1_init,  
  statenames = mod1_s_nm,
  skeleton = vectorfield(f = mod1_skel), 
  params = parms_mod1, 
  paramnames = names(parms_mod1) 
)

# Check that initial conditions sum to 1
stopifnot(sum(rinit(mod1)) == 1)

# Test simulation ----------------------------------------------------------
coef(mod1, names(parms_mod1)) <- unname(parms_mod1)
coef(mod1, "Rzero") <- 15
coef(mod1, c("rho", "alpha")) <- c(2, 1 / 300)
coef(mod1, "pV") <- 0

mod1_tj <- trajectory(object = mod1, format = "data.frame") %>% 
  select(-.id) %>% 
  mutate(N = rowSums(across(.cols = all_of(mod1_s_nm))), 
         lambda = coef(mod1, "Rzero") * (parms_mod1["gamma"] + parms_mod1["mu"]) * I,
         A = 1 / lambda,
         seroprev = Rp1 + Rp2 + Vp, 
         ppv = Rp1 / seroprev)

mod1_tj_long <- mod1_tj %>% 
  pivot_longer(cols = -c("time"), names_to = "var", values_to = "prop")

pl <- ggplot(data = mod1_tj_long %>% filter(time >= 0), 
             mapping = aes(x = time, y = prop)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y") + 
  labs(x = "Time (years)", y = "Proportion")
print(pl)

# Print summary
sumry <- mod1_tj_long %>% 
  filter(time >= max(time) - 19) %>% 
  group_by(var) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup()
print(sumry)

# Calibrate alpha (waning rate) to reach target seroprevalence -----------------------------
coef(mod1, names(parms_mod1)) <- unname(parms_mod1) # Reset parameters
target_prev <- 0.2 # Target seroprevalence
print(coef(mod1))

# Generate list of alpha parameters for given R0 and different rho
parms_df <- tibble(rho = c(0.5, 1), 
                   Rzero = 15,
                   alpha = list(1 / seq(20, 500, by = 10))) %>% 
  unnest(cols = "alpha")

parms_mat <- parmat(params = coef(mod1), nrep = nrow(parms_df))
parms_mat[colnames(parms_df), ] <- t(parms_df)

parms_df <- parms_df %>% 
  mutate(.id = seq_len(nrow(.)), 
         D = 1 / alpha)

# Simulate model
sims_mod1 <- trajectory(object = mod1, 
                        params = parms_mat, 
                        format = "data.frame") 
sims_mod1 <- sims_mod1 %>% 
  mutate(.id = as.integer(.id), 
         seroprev = Rp1 + Rp2 + Vp, 
         ppv = Rp1 / seroprev) %>% 
  pivot_longer(cols = -c(".id", "time"), names_to = "var", values_to = "prop") %>% 
  filter(time >= max(time) - 19, var %in% c("seroprev", "ppv")) %>%
  group_by(.id, var) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup() %>% 
  left_join(parms_df)

# Plot
pl <- ggplot(data = sims_mod1 %>% filter(var == "seroprev"), 
             mapping = aes(x = D, y = 100 * prop, color = factor(rho))) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = 1e2 * target_prev, linetype = "dashed") + 
  theme_classic() +
  labs(x = "Average duration of immunity (years)", 
       y = "Seroprevalence (%)", 
       color = expression(rho), 
       title = "Estimation of alpha, base model")
print(pl)

# Parameter estimates
est_mod1 <- sims_mod1 %>% 
  filter(var == "seroprev") %>% 
  mutate(prop_diff = abs(prop - target_prev)) %>% 
  group_by(rho) %>% 
  mutate(is_best = prop_diff == min(prop_diff)) %>% 
  ungroup() %>% 
  filter(is_best)

print(est_mod1)

tmp <- sims_mod1 %>% 
  filter(.id %in% est_mod1$.id) %>% 
  pivot_wider(values_from = "prop", names_from = "var") %>% 
  mutate(ppv = 1e2 * ppv)
print(tmp %>% mutate(ppv = round(ppv, 0)))

# Extended model: Two V compartments, two R compartments, two boosting coefficients ---------------------------------------------

# Names of state variables
mod2_s_nm <- c("S", "E", "I", "R1", "R2", "Rp1", "Rp2", "V1", "V2", "Vp")

# Deterministic skeleton (ODEs)
mod2_skel <- Csnippet("
    // Force of infection
    double lambda = Rzero * (gamma + mu) * I; // Force of infection
    double v_rate = pV / (1.0 - pV) * mu; // Effective vaccine coverage
    
    // ODEs
    DS = mu + 2.0 * alpha * (R2 + V2) - (lambda + v_rate + mu) * S;
    DE = lambda * S - (sigma + mu) * E;
    DI = sigma * E - (gamma + mu) * I;
    DRp1 = gamma * I - (1.0 / tn + mu) * Rp1;
    DR1 = Rp1 / tn + Rp2 / tn - (rho * lambda + 2.0 * alpha + mu) * R1;
    DR2 = 2.0 * alpha * (R1 - R2) - (kappa * rho * lambda + mu) * R2;
    DRp2 = rho * lambda * (R1 + kappa * R2) - (1.0 / tn + mu) * Rp2;
    DV1 = v_rate * S + Vp / tn - (rho * lambda + 2.0 * alpha + mu) * V1;
    DV2 = 2.0 * alpha  * (V1 - V2) - (kappa * rho * lambda + mu) * V2;
    DVp = rho * lambda * (V1 + kappa * V2) - (1.0 / tn + mu) * Vp;
  ")

# Initializer (based on equilibria of VSIRS model)
mod2_init <- Csnippet("
  S = 1.0 / Rzero; 
  V1 = pV / (1.0 - pV) * mu / (alpha + mu) / Rzero; 
  V2 = 0.0;
  E = 0.0; 
  I = (1.0 - S - V1) / (1.0 + gamma / (alpha + mu));
  R1 = 1.0 - V1 - S - I;
  R2 = 0.0;
  Rp1 = 0.0; 
  Rp2 = 0.0; 
  Vp = 0.0; 
")

# Base model parameters
parms_mod2 <- c("mu" = 1 / 80, # Birth/death rate
                "pV" = 0.9, # Effective vaccination coverage
                "Rzero" = 15, # Basic reproduction number
                "sigma" = 365 / 8, # 1 / average latent period
                "gamma" = 365 / 15, # 1 / average infectious period
                "tn" = 1.9, # Average duration of seropositivity
                "alpha" = 1 / 40, # 1 / average duration of immunity
                "rho" = 0.5, # Boosting coefficient in V1/R1 stage
                "kappa" = 1) # Relative boosting coefficient in V2/R2 stage

# Define the model
mod2 <- pomp(
  data = data.frame(time = seq(1, 5e2, by = 1), cases = NA),
  obsnames = "cases", 
  times = "time",
  t0 = 0, 
  rinit = mod2_init,
  statenames = mod2_s_nm,
  skeleton = vectorfield(f = mod2_skel), 
  params = parms_mod2, 
  paramnames = names(parms_mod2) 
)

# Check that initial conditions sum to 1
stopifnot(sum(rinit(mod2)) == 1)

# Test simulation ----------------------------------------------------------
coef(mod2, names(parms_mod2)) <- unname(parms_mod2)
coef(mod2, "Rzero") <- 15
coef(mod2, c("rho", "alpha", "kappa")) <- c(2, 1 / 40, 1)
# coef(mod2, "pV") <- 0

mod2_tj <- trajectory(object = mod2, format = "data.frame") %>% 
  select(-.id) %>% 
  mutate(N = rowSums(across(.cols = all_of(mod2_s_nm))), 
         lambda = coef(mod2, "Rzero") * (parms_mod2["gamma"] + parms_mod2["mu"]) * I,
         A = 1 / lambda,
         seroprev = Rp1 + Rp2 + Vp, 
         ppv = Rp1 / seroprev)

mod2_tj_long <- mod2_tj %>% 
  pivot_longer(cols = -c("time"), names_to = "var", values_to = "prop")

pl <- ggplot(data = mod2_tj_long %>% filter(time >= 0), 
             mapping = aes(x = time, y = prop)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y") + 
  labs(x = "Time (years)", y = "Proportion")
print(pl)

# Print summary
sumry <- mod2_tj_long %>% 
  filter(time >= max(time) - 19) %>% 
  group_by(var) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup()
print(sumry)

# Calibrate alpha to reach target seroprevalence -----------------------------
coef(mod2, names(parms_mod2)) <- unname(parms_mod2) # Reset parameters
print(coef(mod2))
target_prev <- 0.2 # Target seroprevalence

# List of candidates alphas
parms_df <- tibble(rho = c(0.5, 1), 
                   Rzero = 15,
                   kappa = 1,
                   alpha = list(1 / seq(5, 200, by = 5), 
                                1 / seq(5, 500, by = 5))) %>% 
  unnest(cols = "alpha")

parms_mat <- parmat(params = coef(mod2), nrep = nrow(parms_df))
parms_mat[colnames(parms_df), ] <- t(parms_df)

parms_df <- parms_df %>% 
  mutate(.id = seq_len(nrow(.)), 
         D = 1 / alpha)

# Simulate model
sims_mod2 <- trajectory(object = mod2, 
                        params = parms_mat, 
                        format = "data.frame") 
sims_mod2 <- sims_mod2 %>% 
  mutate(.id = as.integer(.id), 
         seroprev = Rp1 + Rp2 + Vp, 
         ppv = Rp1 / seroprev) %>% 
  pivot_longer(cols = -c(".id", "time"), names_to = "var", values_to = "prop") %>% 
  filter(time >= max(time) - 19, var %in% c("seroprev", "ppv")) %>%
  group_by(.id, var) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup() %>% 
  left_join(parms_df)

# Find best values
est_mod2 <- sims_mod2 %>% 
  filter(var == "seroprev") %>% 
  mutate(prop_diff = abs(prop - target_prev)) %>% 
  group_by(rho) %>% 
  mutate(is_best = prop_diff == min(prop_diff)) %>% 
  ungroup() %>% 
  filter(is_best)
print(est_mod2)

# Plot 
pl <- ggplot(data = sims_mod2 %>% filter(var == "seroprev"), 
             mapping = aes(x = D, y = 1e2 * prop, color = factor(rho))) + 
  geom_line() + 
  geom_point() + 
  geom_hline(yintercept = 1e2 * target_prev, linetype = "dashed") + 
  theme_classic() + 
  labs(x = "Average duration of immunity (years)", 
       y = "Seroprevalence (%)", 
       color = expression(rho), 
       title = "Estimation of alpha, extended model")
print(pl)

tmp <- sims_mod2 %>% 
  filter(.id %in% est_mod2$.id) %>% 
  pivot_wider(values_from = "prop", names_from = "var") %>% 
  mutate(ppv = 1e2 * ppv)
print(tmp %>% mutate(ppv = round(ppv, 0)))

#######################################################################################################
# END
#######################################################################################################
