#######################################################################################################
# Cluster contact matrices
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
source("f-Project_M.R")
source("f-Upsize_M.R")
source("f-PlotMatrix.R")
source("f-compute_R0.R")
library(cluster)
library(factoextra)
library(NbClust)
library(adespatial)
library(clValid)
library(mgcv)
library(gratia)

theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)

# Load country-level social data from Mistry, 2021 ------------------------
path_dat <- "_data/_contact_matrices/Mistry_2021/"
l_files <- list.files(path = path_dat)

# Extract names of countries
country_nm <- map_chr(.x = l_files, .f = \(x) str_extract(string = x, pattern = "^(.*?)(?=_country_level)"))
dat_SCM <- vector(mode = "list", length = length(l_files))
names(dat_SCM) <- country_nm
dat_pop <- matrix(data = 0, nrow = 80, ncol = length(country_nm), 
                  dimnames = list(as.character(0:79), country_nm))

for(i in seq_along(l_files)) {
  file_cur <- l_files[i] # Current file
  dat_SCM[[i]] <- read_csv(file = paste0(path_dat, file_cur), col_names = F, show_col_types = F) %>% 
    as.matrix()
  colnames(dat_SCM[[i]]) <- NULL
  
  # Remove age groups 80-84 yr
  dat_SCM[[i]] <- dat_SCM[[i]][-c(81:85), -c(81:85)]
  
  # Extract empirical data 
  demog_dat <- read_csv(file = sprintf("_data/_demog/_2010/%s_country_level_age_distribution_85.csv", country_nm[i]), 
                        col_names = c("age", "pop"), 
                        col_types = "d") %>% 
    arrange(age) %>% 
    filter(age <= 79)
  stopifnot(nrow(demog_dat) == 80)
  
  # Project to model population structure 
  N_vec <- demog_dat$pop
  dat_SCM[[i]] <- Project_M(M = dat_SCM[[i]], N_tar = N_vec)
  dat_pop[, country_nm[i]] <- N_vec
}

PlotMatrix(M_in = dat_SCM[["France"]], plot_title = "")

# Next-generation matrices ------------------------------------------------
dat_NGM <- vector(mode = "list", length = length(l_files))
names(dat_NGM) <- country_nm

# Transmission parameters
q1_val <- 0.09373881
q2_val <- 0.5117505 * q1_val
q3_val <- 0.1725993 * q2_val
theta_val <- 0.9874919
gamma_val <- 365 / 15
q_vec <- c(rep(q1_val, 10), rep(q2_val, 10), rep(q3_val, 60))

for(i in seq_along(dat_NGM)) {
  nm_cur <- names(dat_NGM)[i]
  
  # Birth rate
  b_rate <- read_csv2(file = "_data/_demog/birth_rates_2010.csv", col_names = T, col_types = "cd")
  b_rate <- b_rate$birth_rate[b_rate$country == nm_cur] / 1e3
  stopifnot(length(b_rate) == 1 && is.numeric(b_rate))
  
  # Aging rates
  delta_vec <- numeric(80)
  N_vec <- dat_pop[, nm_cur]
  Ntot_val <- sum(N_vec)
  delta_vec[1] <- b_rate * Ntot_val / N_vec[1]
  delta_vec[-1] <- delta_vec[1] * N_vec[1] / N_vec[-1]
  
  print(sprintf("Country: %s, Ntot = %.1fM, b_rate=%.3f per yr", nm_cur, Ntot_val / 1e6, b_rate))
  
  # NGM 
  dat_NGM[[nm_cur]] <- compute_R0(theta = theta_val, 
                                  gamma = gamma_val, 
                                  N = N_vec, 
                                  q = q_vec, 
                                  delta = delta_vec, 
                                  Cmat = 365 * dat_SCM[[nm_cur]], 
                                  type = "SIR")
}

R0_vec <- map_dbl(.x = dat_NGM, .f = \(x) x %>% eigen() %>% pluck("values") %>% abs() %>% max())

# Plot demographic structure ----------------------------------------------
dat_pop_long <- reshape2::melt(dat_pop, value.name = "pop")
dat_pop_long <- dat_pop_long %>% 
  rename("country" = "Var2", "age" = "Var1") %>% 
  select(country, age, pop) %>% 
  group_by(country) %>% 
  mutate(pop_frac = pop / sum(pop)) %>% 
  ungroup()

pl <- ggplot(data = dat_pop_long, 
             mapping = aes(x = age, y = pop_frac, group = country)) + 
  geom_line(color = "grey", alpha = 1) + 
  labs(x = "Age (years)", y = "Relative population")
print(pl)

# Plot degree distribution ------------------------------------------------
c_names <- read.table(file = "_data/list_countries.txt") %>% pluck(1)
dat_SCM_sub <- dat_SCM[c_names]

n_ages <- ncol(dat_SCM[[1]]) # No of age groups

deg_dist <- dat_SCM_sub %>% 
  map(.f = \(x) data.frame(age = 0:(n_ages - 1), r_cont = rowSums(x))) %>% 
  bind_rows(.id = "country")

pl <- ggplot(data = deg_dist %>% filter(age <= 85), 
             mapping = aes(x = age, y = r_cont, color = country)) + 
  geom_line() + 
  scale_color_brewer(palette = "Set3") + 
  theme_classic() + 
  labs(x = "Age (years)", y = "Total contact rate (per day)", color = "Country")
print(pl)

r_Newman <- function(M_in) {
  
  # Use notations from https://link.springer.com/article/10.1007/s40314-013-0017-7 
  
  # Renormalize contact matrix 
  M_norm <- M_in / sum(M_in)
  
  a_i <- rowSums(M_norm)
  b_j <- colSums(M_norm)
  
  out <- (sum(diag(M_norm)) - sum(a_i * b_j)) / (1 - sum(a_i * b_j))
  return(out)
}

# Summary statistics
deg_dist_sum <- deg_dist %>% 
  group_by(country) %>% 
  summarise(age_mode = age[which.max(r_cont)], 
            age_mean = sum(age * r_cont) / sum(r_cont), 
            age_var = sum((age ^ 2) * r_cont) / sum(r_cont) - (age_mean ^ 2), 
            age_sd = sqrt(age_var),
            r_cont_tot = sum(r_cont) / 2, 
            mean_r_cont = mean(r_cont), 
            sd_r_cont = sd(r_cont), 
            max_r_cont = max(r_cont)
  ) %>% 
  ungroup() %>% 
  mutate(assort_r = map_dbl(.x = dat_SCM_sub, .f = \(x) r_Newman(x)))

# Estimate age with peak contacts using a GAM -----------------------------
c_name <- "United_States"

dd <- deg_dist %>% 
  filter(country == c_name, age <= 30) %>% 
  arrange(age)

# Fit GAM
mod <- gam(formula = r_cont ~ s(age), data = dd)

# Extract samples of fitted values
fs <- fitted_samples(model = mod, n = 500) %>% 
  mutate(age = .row - 1)

# Plot
pl <- ggplot(data = dd, mapping = aes(x = age, y = r_cont)) + 
  geom_line(data = fs, 
            mapping = aes(x = age, y = .fitted, group = .draw), color = "grey", alpha = 0.5) + 
  geom_line() 
print(pl)

# Estimate peak age
age_max <- fs %>% 
  group_by(.draw) %>% 
  summarise(age_max = age[which.max(.fitted)]) %>% 
  ungroup()
print(quantile(age_max$age_max, probs = c(0.025, 0.5, 0.975)))

# Put data in matrix format for clustering (row: country) --------------------------------
dat_cur <- dat_NGM

dat_cur_mat <- dat_cur %>% 
  map(.f = \(x) as.numeric(x)) %>% 
  bind_rows() %>% 
  t()

# Determine optimal no of clusters ---------------------------------------------

# Visualize dissimilarity matrix
dist_nm <- "manhattan" # Name of distance measure
dist_mat <- get_dist(x = dat_cur_mat, method = dist_nm, stand = F)
fviz_dist(dist.obj = dist_mat, gradient = list(high = "#08306b", mid = "#6baed6", low = "#f7fbff"))

# Tried to run this, but it took too long 
# nb <- NbClust(data = dat_cur_mat,
#               diss = dist_mat,
#               distance = NULL,
#               index = 'silhouette', 
#               min.nc = 5,
#               max.nc = 15,
#               method = "ward.D2")

# For the silhouette method, both pam and hclust suggest an optimal no of 10
nb <- clValid(obj = dat_cur_mat, 
              metric = dist_nm,
              nClust = 5:16, 
              clMethods = c("hierarchical", "agnes", "diana", "fanny", "clara", "pam"), 
              validation = "internal", 
              method = "ward")
print(summary(nb))

# Other function to determine optimal no of clusters
fviz_nbclust(x = dat_cur_mat, 
             #FUNcluster = pam, 
             FUNcluster = hcut, 
             diss = dist_mat, 
             #method = "wss", 
             method = "silhouette", 
             k.max = 20, 
             nboot = 50)

# Run clustering for k = 10 -----------------------------------------------

# Hierarchical clustering
cl_hc <- hclust(d = dist_mat, method = "ward.D2")
cl_agnes <- cluster::agnes(x = dist_mat, diss = T, metric = dist_nm, stand = F, method = "ward")
cl_which <- cl_agnes
k_which <- 12

fviz_dend(x = cl_which, k = k_which, type = "phylogenic", k_colors = "npg")
fviz_dend(x = cl_which, k = k_which, type = "rectangle", k_colors = "npg")

# PAM algorithm
#cl_pam <- pam(x = dist_mat, k = 5, diss = T)

# Cluster by age with contiguity constraints ------------------------------
# dat_age <- dat_SCM[["United-States"]]
# dat_age <- dat_age[26:85, 26:85]
# rownames(dat_age) <- paste0("Age ", 25:84)
# #rownames(dat_age) <- paste0("Age ", 0:84)
# d_age <- get_dist(x = dat_age, method = "manhattan", stand = F)
# fviz_dist(dist.obj = d_age)
# 
# nb <- clValid(obj = dat_age, 
#               metric = "manhattan",
#               nClust = 2:10, 
#               clMethods = c("hierarchical"), 
#               validation = "internal", 
#               method = "ward")
# print(summary(nb))
# 
# cl_age <- constr.hclust(d = d_age, method = "ward.D2", chron = T)
# cutree(tree = cl_age, k = 2)
# plot(cl_age, k = 2)


#######################################################################################################
# END
#######################################################################################################
