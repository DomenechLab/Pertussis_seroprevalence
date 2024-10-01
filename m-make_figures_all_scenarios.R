#######################################################################################################
# Make main figures across all scenarios
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)

# Load seroprevalence estimates from Berbers ------------------------------
library(lme4)

dat_sero <- read_xlsx(path = "_data/seroprevalence_data_Berbers_2021.xlsx")
dat_sero <- dat_sero %>% 
  mutate(country = factor(country), 
         age = factor(age), 
         n_neg = n_tot - n_pos, 
         p_pos = n_pos / n_tot)

# Aggregate data over 40-49 and 50-59
dat_sero_agg <- dat_sero %>% 
  group_by(country) %>% 
  summarise(n_pos = sum(n_pos), 
            n_tot = sum(n_tot), 
            n_neg = n_tot - n_pos,
            p_pos = n_pos / n_tot) %>% 
  ungroup()

# Plot 
pl <- ggplot(data = dat_sero_agg, mapping = aes(x = 1, y = 100 * p_pos)) + 
  geom_boxplot(color = "grey") + 
  geom_point(color = "black") + 
  theme_classic() + 
  labs(x = "", y = "Seroprevalence (%)")
print(pl)

# Fit mixed-effects binomial model
# mod_bin <- glmer(cbind(n_pos, n_neg) ~ age + (1 | country), 
#                  data = dat_sero, 
#                  family = binomial(link = "logit"))


# Load simulations ------------------------------------------------------
l_files <- list.files(path = "_saved", pattern = "DV", include.dirs = T)

mod_preds <- vector(mode = "list", length = length(l_files))
for(i in seq_along(mod_preds)) {
  mod_preds[[i]] <- readRDS(file = paste0("_saved/", l_files[i], "/all-", l_files[i], ".rds"))
}

names(mod_preds) <- l_files %>% str_extract(pattern = "DV_[0-9]+") %>% str_remove(pattern = "DV_")

mod_preds <- mod_preds %>% 
  bind_rows(.id = "DV") %>% 
  filter(age_cat %in% c("[40-50)", "[50,60)"), 
         var_nm %in% c("Rp1", "Rp2", "Vp")) %>% 
  group_by(DV, country, .id, time, var_nm) %>% 
  summarise(n = sum(n), 
            pop = sum(pop)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "var_nm", values_from = "n") %>% 
  mutate(sero_prev = (Rp1 + Rp2 + Vp) / pop, 
         sero_PPV = Rp1 / (Rp1 + Rp2 + Vp), 
         DV = as.numeric(DV), 
         alphaV = 1 / DV)

# Calculate median and 95% quantiles
f_list <- vector(mode = 'list', length = 3)
f_list[[1]] <- function(x) quantile(x = x, probs = 0.025)
f_list[[2]] <- function(x) quantile(x = x, probs = 0.5)
f_list[[3]] <- function(x) quantile(x = x, probs = 0.975)
names(f_list) <- c("q_inf", "q_med", "q_sup")

mod_preds_sumry <- mod_preds %>% 
  group_by(DV, alphaV) %>% 
  summarise(across(c(sero_prev, sero_PPV), f_list)) %>% 
  ungroup()


# Plot --------------------------------------------------------------------
pl <- ggplot(data = mod_preds_sumry, 
             mapping = aes(x = 100 * sero_prev_q_med, y = 100 * sero_PPV_q_med, color = factor(100 * alphaV))) + 
  geom_point() + 
  geom_linerange(mapping = aes(xmin = 100 * sero_prev_q_inf, xmax = 100 * sero_prev_q_sup)) + 
  geom_linerange(mapping = aes(ymin = 100 * sero_PPV_q_inf, ymax = 100 * sero_PPV_q_sup)) + 
  geom_boxplot(data = dat_sero_agg, 
               mapping = aes(x = 100 * p_pos, y = 0, color = NULL), 
               width = 1, outlier.colour = "grey") + 
  #scale_color_viridis(discrete = T, option = "mako", direction = -1, end = 0.8, begin = 0.2) + 
  scale_color_brewer(palette = "Reds") + 
  theme_classic() + 
  theme(legend.position = "top") + 
  labs(x = "Seroprevalence (%)", y = "PPV of seropositivity (%)", 
       color = "Waning rate of vaccine-derived immunity (% per yr)")
print(pl)

ggsave(filename = "_figures/test_main.pdf", plot = pl, width = 7, height = 7)
  
#######################################################################################################
# END 
#######################################################################################################  