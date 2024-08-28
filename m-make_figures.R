#######################################################################################################
# Make figures
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
add_rugs <- F # Should rugs be added for every data point

# Parameters of simulations to plot ---------------------------------------
alphaV <- 0.02 # Waning rate of vaccine-derived immunity
rhoV <- 0.25 # Boosting coefficient of vaccine-derived immunity
vacCov <- 0.9 # Effective vaccine coverage
n_doses <- 2 # No of vaccine doses
nm_file <- sprintf("alphaV_%.2f-rhoV_%.2f-vacCov_%.2f-%ddoses", alphaV, rhoV, vacCov, n_doses)
if(!dir.exists(sprintf("_figures/%s", nm_file))) dir.create(sprintf("_figures/%s", nm_file))

# Load simulations -------------------------------------------------------------
l_countries <- list.dirs(path = "_saved", full.names = F, recursive = F)
sims <- vector(mode = "list", length = length(l_countries))
names(sims) <- l_countries

for(i in seq_along(l_countries)) {
  sims[[i]] <- readRDS(file = sprintf("_saved/%s/alphaV_%.2f-rhoV_%.2f-vacCov_%.2f-%ddoses.rds", 
                                      l_countries[i], alphaV, rhoV, vacCov, n_doses)) %>% 
    pluck("merged_ages")
}

sims <- sims %>% 
  bind_rows(.id = "country") %>% 
  mutate(country = fct_recode(country, 
                              "Czechia" = "Czech",
                              "UK" = "United-Kingdom", 
                              "USA" = "United_States"), 
         age_cat = fct_recode(age_cat, 
                              "0-5 mo" = "[0,0.5)", 
                              "6-11 mo" = "[0.5,1)", 
                              "1-4 yo" = "[1,5)", 
                              "5-9 yo" = "[5,10)", 
                              "10-14 yo" = "[10,15)", 
                              "15-19 yo" = "[15,20)", 
                              "20-24 yo" = "[20,25)", 
                              "25-44 yo" = "[25,45)", 
                              "45-64 yo" = "[45,65)", 
                              "65-79 yo" = "[65,Inf]"
         ))

# Subset simulated data -------------------------------------------------------------
sims_cur <- sims %>% 
  filter(time >= max(time) - 19, 
         var_nm %in% c("seroPrev", "seroPPV", "seroInc", "trueInc"), 
         age_cat %in% c("25-44 yo", "45-64 yo", "65-79 yo")) %>% 
  mutate(prop = n / pop)

# Summary: median of all variables
sims_sumry <- sims_cur %>% 
  group_by(country, var_nm, age_cat) %>% 
  summarise(med_prop = median(prop), 
            med_n = median(n)) %>% 
  ungroup()

# Order countries alphabetically
country_order <- sort(levels(sims$country), decreasing = F)

# Recode age groups
sims_cur <- sims_cur %>% 
  mutate(country = factor(x = country, levels = country_order))

# Main figure: seroprevalence (left panel) and PPV (right panel) -------------------------------------------------------------

tmp <- sims_cur %>% 
  filter(var_nm %in% c("seroPrev", "seroPPV")) %>% 
  mutate(val = ifelse(var_nm == "seroPrev", prop, n), 
         var_nm = factor(var_nm, levels = c("seroPrev", "seroPPV"), labels = c("Seroprevalence", "PPV of seropositivity")))

levels(tmp$country) <- rev(levels(tmp$country))

pl <- ggplot(data = tmp, 
             mapping = aes(x = 1e2 * val, y = country, fill = age_cat, color = age_cat)) + 
  geom_density_ridges(scale = 0.9,
                      stat = "density_ridges", 
                      quantile_lines = T, 
                      quantiles = 2, 
                      jittered_points = add_rugs,
                      point_shape = "|", 
                      position = position_points_jitter(height = 0),
                      alpha = 0.5) + 
  facet_wrap(~ var_nm, scales = "free_x", ncol = 2) + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  scale_fill_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "Value (%)", y = "SCM country source", fill = "Age group", color = "") +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        legend.position = "top", 
        strip.text = element_text(size = rel(1))
        ) + 
  guides(color = F)
print(pl)

ggsave(filename = sprintf("_figures/%s/main_Sp_PPV.pdf", nm_file), 
       plot = pl, width = 8, height = 6)

# Sup figure: seroprevalence vs. sero-incidence ---------------------------------------
tmp <- sims_cur %>% 
  select(-c(var_type, prop)) %>% 
  pivot_wider(names_from = "var_nm", values_from = "n")

alpha_val <- 0.5

# Seroprevalence vs. sero-incidence
pl1 <- ggplot(data = tmp, 
              mapping = aes(x = seroPrev / pop, y = seroInc / pop, color = age_cat)) + 
  geom_point(alpha = alpha_val) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "Seroprevalence", y = "Sero-incidence rate (per year)", color = "Age group")
print(pl1)

# PPV vs. true incidence / sero-incidence
pl2 <- ggplot(data = tmp, 
              mapping = aes(x = seroPPV, y = trueInc / seroInc, color = age_cat)) + 
  geom_point(alpha = alpha_val) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "PPV of seropositivity", y = "True incidence rate / sero-incidence rate", color = "Age group")
print(pl2)

pl_all <- pl1 + pl2 + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme_classic() + theme(strip.background = element_blank(), 
                          strip.text = element_text(size = rel(1)),
                          legend.position = "top", 
                          plot.tag = element_text(size = 12)) 
print(pl_all)

ggsave(filename = sprintf("_figures/%s/sup_Seroprev_Seroinc.pdf", nm_file), 
       plot = pl_all, width = 8, height = 6)

# Sup figure: age distribution of cases -----------------------------------

age_dist <- sims %>% 
  filter(var_nm %in% c("Ci1", "Ci2"), 
         time >= max(time) - 19) %>% 
  group_by(country, .id, var_nm, age_cat) %>% 
  summarise(inc_mean = mean(n / pop), 
            inc_sd = sd(n / pop)) %>% 
  ungroup() %>% 
  mutate(var_nm = fct_recode(var_nm, 
                             "Primary infection" = "Ci1", 
                             "Secondary infection" = "Ci2"))

pl <- ggplot(data = age_dist, 
             mapping = aes(x = age_cat, y = 1e5 * inc_mean, color = country, group = interaction(country, .id))) + 
  geom_line(alpha = 1) + 
  scale_color_brewer(palette = "Paired") + 
  facet_wrap(~ var_nm, nrow = 2, scales = "free") + 
  theme_classic() + theme(strip.background = element_blank(), 
                          strip.text = element_text(size = rel(1)),
                          legend.position = c(0.5, 0.8)) + 
  #scale_y_sqrt() + 
  labs(x = "Age group", y = "Incidence rate (per year per 100,000)", color = "")
print(pl)

ggsave(filename = sprintf("_figures/%s/sup_age_distribution.pdf", nm_file), 
       plot = pl, width = 8, height = 8)


#######################################################################################################
# End
#######################################################################################################