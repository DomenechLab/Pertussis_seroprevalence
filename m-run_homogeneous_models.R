#######################################################################################################
# Run simulations of homogeneous SIRWS model with vaccination
# Time unit: years; rates are per year
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
library(pomp)
library(metR)
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)

# Csnippet for ODE models --------------------------------------------------

names_mod <- c("SIRWS", "SIRWS_bis", "SIRRp")
skel <- vector(mode = "list", length = 3)
names(skel) <- names_mod

skel[["SIRWS"]] <- Csnippet("
  double b = R0 * (gamma + mu); 
  double lambda = b * I; 
  
  DV = mu * p_V + rho_V * lambda * W_V - (2.0 * alpha_V + mu) * V; 
  DW_V = 2.0 * alpha_V * V - (rho_V * lambda + 2.0 * alpha_V + mu) * W_V; 
  DS = mu * (1.0 - p_V) + 2.0 * alpha_V * W_V + 2.0 * alpha_R * W_R - (lambda + mu) * S; 
  DI = lambda * S - (gamma + mu) * I; 
  DR = gamma * I + rho_R * lambda * W_R - (2.0 * alpha_R + mu) * R; 
  DW_R = 2.0 * alpha_R * R - (2.0 * alpha_R + rho_R * lambda + mu) * W_R; 
")

skel[["SIRWS_bis"]] <- Csnippet("
  double b = R0 * (gamma + mu); // Transmission rate
  double lambda = b * I; // Force of infection
  double alpha_prime = 1.0 / (0.5 / alpha_R - 1.0 / alpha_P); // Rate of transition R->W_R
  
  DV = mu * p_V + alpha_P * Vp - (rho_V * lambda + alpha_V + mu) * V; 
  DVp = rho_V * lambda * V - (alpha_P + mu) * Vp; 
  DS = mu * (1.0 - p_V) + alpha_V * V + 2.0 * alpha_R * W_R - (lambda + mu) * S; 
  DI = lambda * S - (gamma + mu) * I; 
  DRp = gamma * I + rho_R * lambda * W_R - (alpha_P + mu) * Rp;
  DR = alpha_P * Rp - (alpha_prime + mu) * R; 
  DW_R = alpha_prime * R - (2.0 * alpha_R + rho_R * lambda + mu) * W_R; 
")

skel[["SIRRp"]] <- Csnippet("
  double b = R0 * (gamma + mu); 
  double lambda = b * I; 
  
  DV = mu * p_V + alpha_P * Vp - (rho_V * lambda + alpha_V + mu) * V; 
  DVp = rho_V * lambda * V - (alpha_P + mu) * Vp; 
  DS = mu * (1.0 - p_V) + alpha_V * V + alpha_R * R - (lambda + mu) * S; 
  DI = lambda * S - (gamma + mu) * I; 
  DRp1 = gamma * I - (alpha_P + mu) * Rp1;
  DR = alpha_P * Rp1 + alpha_P * Rp2 - (rho_R * lambda + alpha_R + mu) * R; 
  DRp2 = rho_R * lambda * R - (alpha_P + mu) * Rp2; 
")

# Parameters -------------------------------------------------------
theta <- vector(mode = "list", length = 3)
names(theta) <- names_mod

theta[["SIRWS"]] <-  c(
  "mu" = 1 / 75, # Birth/death rate 
  "R0" = 16, # Basic reproduction number
  "gamma" = 1 / (3.7 / 52), # Recovery rate
  "p_V" = 0, # Effective vaccine coverage
  "alpha_V" = 1 / 34, # Rate of waning vaccine immunity (without boosting)
  "rho_V" = 0.25, # Boosting coefficient of vaccine immunity
  "alpha_R" = 1 / 34, # Rate of waning infection-derived immunity (without boosting)
  "rho_R" = 6.6 # Boosting coefficient of infection-derived immunity
)

theta[["SIRRp"]] <- theta[["SIRWS_bis"]] <- theta[["SIRWS"]]
theta[["SIRRp"]] <- c(theta[["SIRRp"]], "alpha_P" = 1 / 1)
theta[["SIRWS_bis"]] <- c(theta[["SIRWS_bis"]], "alpha_P" = 1 / 1)

# Names of state variables
state_nm <- vector(mode = "list", length = 3)
names(state_nm) <- names_mod

state_nm[["SIRWS"]] <-  c("V", "W_V", "S", "I", "R", "W_R")
state_nm[["SIRWS_bis"]] <- c("V", "Vp", "S", "I", "Rp", "R", "W_R")
state_nm[["SIRRp"]] <- c("V", "Vp", "S", "I", "Rp1", "R", "Rp2")

for(i in seq_along(theta)) {
  theta[[i]] <- c(theta[[i]], setNames(object = rep(0, length(state_nm[[i]])), 
                                       nm = paste0(state_nm[[i]], ".0")))
  theta[[i]]["V.0"] <- unname(theta[["SIRWS"]]["p_V"])
  theta[[i]]["S.0"] <- unname(1 - theta[["SIRWS"]]["p_V"] - 1e-3)
  theta[[i]]["I.0"] <- 1e-3
}

# Transmission rate
beta_val <- unname(theta[[1]]["R0"] * (theta[[1]]["gamma"] + theta[[1]]["mu"]))

# Create POMP models ------------------------------------------------------
mod <- vector(mode = "list", length = 3)
names(mod) <- names_mod

for(i in seq_along(mod)) {
  
  mod[[i]] <- pomp(data = data.frame(time = seq(0, 500, by = 0.05), X = NA),
                   times = "time",
                   t0 = 0,
                   obsnames = "X",
                   statenames = state_nm[[i]],
                   paramnames = names(theta[[i]]), 
                   params = theta[[i]],
                   skeleton = vectorfield(f = skel[[i]]))
  stopifnot(sum(rinit(mod[[i]])) == 1)
  print(i)
}

# Test simulation ----------------------------------------------------------------
for(i in seq_along(mod)) {
  
  tj <- trajectory(object = mod[[i]], format = "data.frame") 
  tj$N <- rowSums(tj[, state_nm[[i]]])
  
  tj <- tj %>% 
    mutate(lambda = beta_val * I)
  
  tj_long <- tj %>% 
    pivot_longer(cols = -c(".id", "time"), names_to = "var", values_to = "prop")
  
  pl <- ggplot(data = tj_long, mapping = aes(x = time, y = prop)) + 
    geom_line() + 
    facet_wrap(~ var, scales = "free_y") +
    scale_y_sqrt() + 
    labs(x = "Time (years)", y = "Proportion", title = names(mod)[i])
  print(pl)
  
  print(filter(tj_long, time == max(time)))
}

# Simulate SIRRp model ----------------------------------------------------
p_V <- 0.9
which_mod <- "SIRRp"
coef(mod[[which_mod]], names(theta[[which_mod]])) <- unname(theta[[which_mod]])
coef(mod[[which_mod]], "p_V") <- p_V
coef(mod[[which_mod]], "alpha_P") <- 1 / 1
coef(mod[[which_mod]], c("V.0", "S.0", "I.0")) <- c(p_V, 1 - p_V - 1e-3, 1e-3)
coef(mod[[which_mod]], c("rho_R", "rho_V")) <- 2
coef(mod[[which_mod]], c("alpha_R", "alpha_V")) <- 1 / 40
p_vec <- coef(mod[[which_mod]])

stopifnot(sum(rinit(mod[[which_mod]])) == 1)

tj <- trajectory(object = mod[[which_mod]], format = "data.frame") 
tj$N <- rowSums(tj[, state_nm[[which_mod]]])

tj <- tj %>% 
  mutate(lambda = beta_val * I, 
         Sp = Vp + Rp1 + Rp2)

tj_long <- tj %>% 
  pivot_longer(cols = -c(".id", "time"), names_to = "var", values_to = "prop")

pl <- ggplot(data = tj_long, mapping = aes(x = time, y = prop)) + 
  geom_line() + 
  facet_wrap(~ var, scales = "free_y") +
  scale_y_sqrt() + 
  labs(x = "Time (years)", y = "Proportion", title = names(mod[[i]]))
print(pl)

print(filter(tj_long, time == max(time)))
lambda_eq <- tj$lambda[which.max(tj$time)]

f_D <- function(rho, alpha, alpha_P, lambda) (alpha_P + rho * lambda) * (alpha + rho * lambda) / (alpha_P * alpha^2)
D <- f_D(rho = p_vec["rho_V"], alpha = p_vec["alpha_V"], alpha_P = p_vec["alpha_P"], lambda = lambda_eq)
print(lambda_eq)
print(D)

# Grid search to calibrate rho_R -------------------------------------------------------------
all_parms <- expand_grid(rho_R = seq(5, 20, length.out = 20), 
                         alpha_R = 1 / seq(10, 100, length.out = 20)) %>% 
  mutate(.id = seq_len(nrow(.)))

p_mat <- parmat(params = coef(mod[["SIRRp"]]), 
                nrep = nrow(all_parms))
p_mat["rho_R", ] <- all_parms$rho_R
p_mat["alpha_R", ] <- all_parms$alpha_R

tj_multi <- bake(file = "_saved/other/simulations_SIRRp_model.rds", 
                 expr = {
                   trajectory(object = mod[["SIRRp"]], 
                              params = p_mat,
                              format = "data.frame") %>% 
                     mutate(.id = as.integer(.id), 
                            Sp = Rp1 + Rp2,
                            lambda = beta_val * I) %>% 
                     left_join(y = all_parms)
                 })
  

tj_multi_long <- tj_multi %>% 
  pivot_longer(cols = -c(".id", "time", "rho_R", "alpha_R"), names_to = "var", values_to = "prop")

eq_vals <- tj_multi %>% filter(time == max(time))

pl <- ggplot(data = eq_vals, mapping = aes(x = rho_R, y = 1 / alpha_R, z = lambda - 0.231)) + 
  geom_tile(mapping = aes(fill = lambda - 0.231)) + 
  geom_contour(color = "black") + 
  geom_text_contour(color = "black") +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") + 
  #scale_fill_scico(palette = "vik", midpoint = 0) + 
  labs(x = expression(rho[R]), y = expression(D[R]))
print(pl)

pl <- ggplot(data = eq_vals, mapping = aes(x = rho_R, y = 1 / alpha_R, z = 100 * Sp)) + 
  geom_tile(mapping = aes(fill = 100 * Sp)) + 
  geom_contour(color = "white") + 
  geom_text_contour(color = "white") +
  scale_fill_viridis(option = "viridis") + 
  labs(x = expression(rho[R]), y = expression(D[R]))
print(pl)

# Estimate rho_R and alpha_R to reach target FoI and seroprevalence ----------------------------------------------
coef(mod[["SIRRp"]], names(theta[["SIRRp"]])) <- unname(theta[["SIRRp"]])

dist_fun <- function(par, FoI_tar, Sp_tar, cut_off = 125) {
  # Args: 
  # par: parameters to estimate (rho_R, alpha_R) [2-vector]
  # FoI_tar: target force of infection (per year) [scalar]
  # Sp_tar: target seroprevalence [scalar]
  # Returns: squared distance between simulated and target FoI and Sp
  
  # Set parameters
  rho_R_test <- exp(par[1])
  D_R_test <- exp(par[2])
  
  if(cut_off == 125) {
    coef(mod[["SIRRp"]], "alpha_P") <- 1 / 0.75
  } else if(cut_off == 50) {
    coef(mod[["SIRRp"]], "alpha_P") <- 1 / 2.42
  }
  
  coef(mod[["SIRRp"]], c("alpha_R", "rho_R")) <- c(1 / D_R_test, rho_R_test)
  
  # Run model until equilibrium
  sim_cur <- trajectory(object = mod[["SIRRp"]], format = "data.frame") %>% 
    mutate(Sp = Rp1 + Rp2, 
           FoI = beta_val * I)
  
  # Extract equilibrium
  Sp_pred <- sim_cur$Sp[which.max(sim_cur$time)]
  FoI_pred <- sim_cur$FoI[which.max(sim_cur$time)]
  
  out <- ((Sp_pred - Sp_tar) ^ 2) + ((FoI_pred - FoI_tar) ^ 2)
  return(out)
}

x <- dist_fun(par = log(c(5, 10)), FoI_tar = 0.231, Sp_tar = 0.5, cut_off = 125)

fit <- optim(par = log(c(5, 10)), 
             fn = dist_fun, 
             FoI_tar = 0.257, 
             Sp_tar = 0.238, 
             cut_off = 125,
             method = "Nelder-Mead")
stopifnot(fit$convergence == 0)

par_est <- exp(fit$par)
print(par_est %>% round(1))
print(fit$value)

#######################################################################################################
# END
#######################################################################################################