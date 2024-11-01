#######################################################################################################
# Make main figures across all scenarios
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
library(lemon)
library(latex2exp)
library(ggrepel)
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
t_V <- 150 # Year of vaccine introduction
path_to_res <- "_saved/estimation_6_countries/rho_5.00/" # Path to saved results
suffix <- str_extract(string = path_to_res, pattern = "rho_[0-9].[0-9]+") # String to append to PDF file name

# Load seroprevalence data ------------------------------------------------
dat_seroprev <- read_xlsx(path = "_data/_seroprevalence/seroprevalence_data_Pebody_Wehlin.xlsx", 
                          sheet = "data", 
                          col_names = T) %>% 
  left_join(y = data.frame(country = c("Finland", "France", "Germany", "Netherlands", "Lithuania", "Romania"), 
                           country_short = c("FI", "FR", "DE", "NL", "LT", "RO"))) %>% 
  select(country, country_short, everything()) %>% 
  mutate(age_cat = factor(age_cat), 
         age_cat = fct_recode(age_cat, "[65,Inf]" = "65+"), 
         binom_test = map2(.x = n_pos, .y = n_tot, .f = \(x, y) binom.test(x = x, n = y, conf.level = 0.95))) %>% 
  mutate(Sp_obs_est = map_dbl(.x = binom_test, .f = ~ pluck(.x, "estimate")), 
         Sp_obs_se = sqrt(Sp_obs_est * (1 - Sp_obs_est) / n_tot),
         Sp_obs_inf = map_dbl(.x = binom_test, .f = ~ pluck(.x, "conf.int", 1)), 
         Sp_obs_sup = map_dbl(.x = binom_test, .f = ~ pluck(.x, "conf.int", 2)))

# Load simulated data -----------------------------------------------------

# Load simulation details for every country
sims_details <- read_xlsx(path = "_data/_seroprevalence/simulation_details.xlsx", sheet = "data", col_names = T)


# List directories (each directory corresponds to a different rate of vaccine-derived immunity)
l_dirs <- list.files(path = path_to_res, pattern = "DV", include.dirs = T)
print(l_dirs)

# Store results into a named list (name: average duration of vaccine-derived immunity)
mod_preds <- vector(mode = "list", length = length(l_dirs))
names(mod_preds) <- l_dirs %>% str_extract(pattern = "DV_[0-9]+") %>% str_remove(pattern = "DV_")

for(i in seq_along(mod_preds)) {
  
  dir_cur <- l_dirs[i]
  
  # List files in current directory
  l_files <- list.files(path = paste0(path_to_res, dir_cur), pattern = ".rds")
  country_nm <- str_extract(string = l_files, pattern = "^([^-]+)") # Extract country names
  
  mod_preds[[i]] <- vector(mode = "list", length = length(l_files))
  names(mod_preds[[i]]) <- country_nm
  
  # Load data
  for(j in seq_along(l_files)) {
    file_cur <- l_files[j]
    country_cur <- country_nm[j]
    mod_preds[[i]][[j]] <- readRDS(file = sprintf("%s%s/%s", path_to_res, dir_cur, file_cur)) %>% 
      pluck("merged_ages") %>% 
      filter(between(time, 
                     sims_details$t1[sims_details$country == country_cur] + t_V, 
                     sims_details$t2[sims_details$country == country_cur] + t_V), 
             age_cat %in% dat_seroprev$age_cat)
  }
  
  # Transform list into data frame
  mod_preds[[i]] <- mod_preds[[i]] %>% 
    bind_rows(.id = "country")
}

mod_preds <- mod_preds %>% 
  bind_rows(.id = "DV")

# Calculate seroprevalence and PPV ----------------------------------------
mod_preds_sub <- mod_preds %>% 
  filter(var_nm %in% c("Rp1", "Rp2", "Vp")) %>% 
  group_by(DV, country, .id, time, age_cat, var_nm) %>% 
  summarise(n = sum(n), 
            pop = sum(pop)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "var_nm", values_from = "n") %>% 
  mutate(Sp_pred = (Rp1 + Rp2 + Vp) / pop, 
         PPV_pred = Rp1 / (Rp1 + Rp2 + Vp), 
         DV = as.numeric(DV), 
         alphaV = 1 / DV)

# Calculate median and 95% quantiles
f_list <- vector(mode = 'list', length = 3)
f_list[[1]] <- function(x) quantile(x = x, probs = 0.025)
f_list[[2]] <- function(x) quantile(x = x, probs = 0.5)
f_list[[3]] <- function(x) quantile(x = x, probs = 0.975)
names(f_list) <- c("q_inf", "q_med", "q_sup")

mod_preds_sumry <- mod_preds_sub %>% 
  filter(!is.na(PPV_pred)) %>% 
  group_by(DV, alphaV, country, age_cat) %>% 
  summarise(across(c(Sp_pred, PPV_pred), f_list)) %>% 
  ungroup()

# Merge simulated and observed data
pred_obs <- mod_preds_sumry %>% 
  left_join(y = dat_seroprev %>% select(-c(binom_test))) %>% 
  mutate(age_cat_no = ifelse(age_cat %in% c("[20,40)", "[20,45)"), 1, ifelse(age_cat == "[65,Inf]", 3, 2)), 
         DV = as.factor(DV), 
         DV_title = fct_relabel(.f = DV, .fun = ~ paste0("Dv = ", .x, " years"))) %>% 
  select(DV, DV_title, alphaV, country, country_short, age_cat, age_cat_no, everything())

# Estimation performance
est_perf <- pred_obs %>% 
  filter(!is.na(Sp_obs_est), age_cat_no == 1) %>% 
  group_by(DV, DV_title, alphaV, age_cat_no) %>% 
  summarise(rmse = sqrt(mean((Sp_pred_q_med - Sp_obs_est) ^ 2)), # Root mean squared error
            mab = mean(abs(Sp_pred_q_med - Sp_obs_est)), # Mean absolute boas
            marb = mean(abs(Sp_pred_q_med - Sp_obs_est) / Sp_obs_est),  # Mean absolute relative bias
            me_mod = list(lm(formula = Sp_obs_est ~ 1 + offset(Sp_pred_q_med), # Weighted linear model for mean error
                             weights = 1 / (Sp_obs_se ^ 2)))
  ) %>% 
  ungroup() %>% 
  mutate(me_est = map_dbl(.x = me_mod, .f = ~ -summary(.x)$coefficients[1, 1]), # -Estimate of intercept (model - obs)
         me_se = map_dbl(.x = me_mod, .f = ~ summary(.x)$coefficients[1, 2]), # SE of intercept
         me_pval = map_dbl(.x = me_mod, .f = ~ summary(.x)$coefficients[1, 4])) # P-value for null value of intercept

# DV values to include: range with ns difference from 0 OR first value with significant positive/negative value for MWSE
est_perf$include <- (est_perf$alphaV %in% range(est_perf$alphaV[est_perf$me_pval >= 0.05])) | # Min-Max of alphaV values with non-NS MWSE
  (est_perf$alphaV == min(est_perf$alphaV[est_perf$me_est > 0 & est_perf$me_pval < 0.05])) | # Min alphaV value with S positive MWSE 
  (est_perf$alphaV == max(est_perf$alphaV[est_perf$me_est < 0 & est_perf$me_pval < 0.05])) # Max alphaV value with S significant MWSE 

stopifnot(sum(est_perf$include) == 4)
print(sprintf("Best DV=%.0f yr, best MWSE=%.1f %%", 
              1 / est_perf$alphaV[which.min(abs(est_perf$me_est))], 
              100 * est_perf$me_est[which.min(abs(est_perf$me_est))]))

print("PPV range")
tmp <- pred_obs %>% 
  filter(DV %in% est_perf$DV[est_perf$me_pval >= 0.05], age_cat_no == 1)
print(range(100 * tmp$PPV_pred_q_med) %>% round(0))

# Weighted linear model ---------------------------------------------------
dat_test <- pred_obs %>% 
  filter(DV == 10, age_cat_no == 1) %>% 
  arrange(country) %>% 
  mutate(w = 1 / (Sp_obs_se ^ 2), 
         w_norm = w / sum(w))

mod <- lm(formula = Sp_obs_est ~ 1 + offset(Sp_pred_q_med), 
          data = dat_test, 
          weights = w_norm)
#print(summary(mod))

# Plot weights
pl <- ggplot(data = dat_test, mapping = aes(x = country_short, y = w)) + 
  geom_col() + 
  geom_text(color = "blue", mapping = aes(y = w + 5e3, label = paste0("n = ", n_tot))) + 
  labs(x = "Country", y = "Weigth")
print(pl)

# Make plot in first adult age group--------------------------------------------------------------
# See https://sahirbhatnagar.com/blog/2016/02/08/ggplot2-facet-wrap-labels/ for how to add the mathematical annotations
appender <- function(string) TeX(paste0("$\\alpha_R^{-1} = \\alpha_V^{-1}= $", string, " yr"))

pl <- ggplot(data = pred_obs %>% filter(age_cat_no == 1, DV %in% est_perf$DV[est_perf$include]), 
             mapping = aes(x = 100 * Sp_obs_est, 
                           y = 100 * Sp_pred_q_med, 
                           label = country_short)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + 
  geom_linerange(mapping = aes(xmin = 100 * Sp_obs_inf, xmax = 100 * Sp_obs_sup), color = "grey") + 
  geom_linerange(mapping = aes(ymin = 100 * Sp_pred_q_inf, ymax = 100 * Sp_pred_q_sup), color = "grey") + 
  geom_point(mapping = aes(colour = 100 * PPV_pred_q_med, shape = cut_off), size = rel(3)) + 
  scale_color_viridis(discrete = F, option = "magma", direction = 1, begin = 0, end = 1, trans = "identity") + 
  scale_shape_manual(values = c(17, 16)) + 
  geom_text_repel(force_pull = 0) + 
  geom_text(data = est_perf %>% mutate(x = 5.5, y = 0) %>% filter(age_cat_no == 1, include), 
            mapping = aes(x = x, y = y, 
                          label = sprintf("MWSE = %.1f%% (%s)", 100 * me_est, if_else(me_pval < 0.05, "*", "ns"))), 
            size = 11, size.unit = "pt") + 
  facet_rep_wrap(~ as.character(round(1 / alphaV, 0)), 
                 labeller = as_labeller(x = appender, default = label_parsed), 
                 scales = "fixed", 
                 dir = "h", 
                 ncol = 2) + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid = element_blank(),
        #axis.line = element_line(color = "black", size = 0.5),  # Keep axis 
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) + 
  labs(x = "Observed seroprevalence (%)", 
       y = "Predicted seroprevalence (%)", 
       color = "PPV (%)", 
       shape = "IgG seropositivity threshold")
print(pl)

ggsave(filename = sprintf("_figures/main/fig_empirical_comparison-%s.pdf", suffix), 
       plot = pl, width = 8, height = 8)

#######################################################################################################
# END 
#######################################################################################################  