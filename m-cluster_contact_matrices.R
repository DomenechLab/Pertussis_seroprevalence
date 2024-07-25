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

theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)

# Load country-level social data from Mistry, 2021 ------------------------
# path_from <- "_data/_contact_matrices/Mistry_2021/_all/"
# path_to <- "_data/_contact_matrices/Mistry_2021/_country_level/"
# l_files <- list.files(path = path_from, pattern = "country_level_M_overall_contact_matrix_85") 
# 
# for(f in l_files) {
#   file.copy(from = paste0(path_from, f), to = paste0(path_to, f))
# }

path_dat <- "_data/_contact_matrices/Mistry_2021/"
l_files <- list.files(path = path_dat)

# Extract names of countries
country_nm <- map_chr(.x = l_files, .f = \(x) str_split_i(string = x, pattern = "_", i = 1))
country_nm[country_nm == "United"] <- "United-States"
country_nm[country_nm == "South"] <- "South-Africa"

dat_SCM <- vector(mode = "list", length = length(l_files))
names(dat_SCM) <- country_nm

for(i in seq_along(l_files)) {
  file_cur <- l_files[i] # Current file
  dat_SCM[[i]] <- read_csv(file = paste0(path_dat, file_cur), col_names = F) %>% 
    as.matrix()
  colnames(dat_SCM[[i]]) <- NULL
}

PlotMatrix(M_in = dat_SCM[["France"]], plot_title = "")


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
              metric = "manhattan",
              nClust = 5:20, 
              clMethods = c("hierarchical", "pam", "kmeans"), 
              validation = "internal", 
              method = "ward")
print(summary(nb))

# Other function to determine optimal no of clusters
fviz_nbclust(x = dat_SCM_mat, 
             FUNcluster = pam, 
             diss = dist_mat, 
             method = "silhouette", 
             k.max = 20, 
             nboot = 50)

# Run clustering for k = 10 -----------------------------------------------

# Hierarchical clustering
cl_hc <- hclust(d = dist_mat, method = "ward.D2")
fviz_dend(x = cl_hc, k = 10, type = "phylogenic")

# PAM algorithm
cl_pam <- pam(x = dist_mat, k = 10, diss = T)


# Cluster by age with contiguity constraints ------------------------------
dat_age <- dat_SCM[["United-States"]]
dat_age <- dat_age[26:85, 26:85]
rownames(dat_age) <- paste0("Age ", 25:84)
#rownames(dat_age) <- paste0("Age ", 0:84)
d_age <- get_dist(x = dat_age, method = "manhattan", stand = F)
fviz_dist(dist.obj = d_age)

nb <- clValid(obj = dat_age, 
              metric = "manhattan",
              nClust = 2:10, 
              clMethods = c("hierarchical"), 
              validation = "internal", 
              method = "ward")
print(summary(nb))

cl_age <- constr.hclust(d = d_age, method = "ward.D2", chron = T)
cutree(tree = cl_age, k = 2)
plot(cl_age, k = 2)


#######################################################################################################
# END
#######################################################################################################
