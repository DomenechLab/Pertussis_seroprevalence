#######################################################################################################
# Cluster contact matrices
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
source("f-Project_M.R")
source("f-Upsize_M.R")
source("f-PlotMatrix.R")
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
}

PlotMatrix(M_in = dat_SCM[["France"]], plot_title = "")

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
dat_SCM_mat <- dat_SCM %>% 
  map(.f = \(x) as.numeric(x)) %>% 
  bind_rows() %>% 
  t()

# Determine optimal no of clusters ---------------------------------------------

# Visualize dissimilarity matrix
dist_nm <- "manhattan" # Name of distance measure
dist_mat <- get_dist(x = dat_SCM_mat, method = dist_nm, stand = F)
fviz_dist(dist.obj = dist_mat)

# Tried to run this, but it took too long 
# nb <- NbClust(data = dat_SCM_mat, 
#               distance = "canberra", 
#               diss = NULL, 
#               min.nc = 5, 
#               max.nc = 15,  
#               method = "ward.D2")

# For the silhouette method, both pam and hclust suggest an optimal no of 10
nb <- clValid(obj = dat_SCM_mat, 
              metric = dist_nm,
              nClust = 5:20, 
              clMethods = c("hierarchical", "pam"), 
              validation = "internal", 
              method = "ward")
print(summary(nb))

# Other function to determine optimal no of clusters
fviz_nbclust(x = dat_SCM_mat, 
             #FUNcluster = pam, 
             FUNcluster = hcut, 
             diss = dist_mat, 
             #method = "gap_stat", 
             method = "silhouette", 
             k.max = 20, 
             nboot = 50)

# Run clustering for k = 10 -----------------------------------------------

# Hierarchical clustering
cl_hc <- hclust(d = dist_mat, method = "ward.D2")
plot(cl_hc)
fviz_dend(x = cl_hc, k = 10, type = "phylogenic")
fviz_dend(x = cl_hc, k = 10, type = "rectangle", horiz = F)

# PAM algorithm
cl_pam <- pam(x = dist_mat, k = 10, diss = T)

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
