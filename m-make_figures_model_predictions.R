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
DV <- 50 # Average duration of vaccine-derived immunity (in years)
DR <- 50 # Average duration of infection-derived immunity (in years)
rhoV <- 2 # Boosting coefficient of vaccine-derived immunity
rhoR <- 2 # Boosting coefficient of infection-derived immunity
vacCov <- 0.9 # Effective vaccine coverage
n_doses <- 3 # No of vaccine doses
#nm_dir <- sprintf("DV_%.0f-rhoV_%.2f-DR_%.0f-rhoR_%.2f-vacCov_%.2f-%ddoses", DV, rhoV, DR, rhoR, vacCov, n_doses)
nm_dir <- sprintf("DV_%.0f-rhoV_%.2f-DR_%.0f-rhoR_%.2f", DV, rhoV, DR, rhoR)
stopifnot(dir.exists(sprintf("_saved/%s", nm_dir)))
if(!dir.exists(sprintf("_figures/%s", nm_dir))) dir.create(sprintf("_figures/%s", nm_dir))
if(!dir.exists(sprintf("_figures/%s/_all", nm_dir))) dir.create(sprintf("_figures/%s/_all", nm_dir))

# Load simulations -------------------------------------------------------------
l_files <- list.files(path = sprintf("_saved/%s", nm_dir), pattern = ".rds")
l_files <- l_files[str_detect(string = l_files, pattern = "all", negate = T)]
l_countries <- str_extract(string = l_files, pattern = "^(.*?)(?=-DV)")

sims <- vector(mode = "list", length = length(l_files))
names(sims) <- l_countries

for(i in seq_along(l_files)) {
  sims[[i]] <- readRDS(file = sprintf("_saved/%s/%s", nm_dir, l_files[i])) %>% 
    pluck("merged_ages")
  
  # Move PDF files
  nm_pdf <- str_replace(l_files[i], ".rds", ".pdf")
  if(file.exists(sprintf("_saved/%s/%s", nm_dir, nm_pdf))) {
    file.copy(from = sprintf("_saved/%s/%s", nm_dir, nm_pdf), 
              to = sprintf("_figures/%s/_all/%s", nm_dir, nm_pdf))
    file.remove(sprintf("_saved/%s/%s", nm_dir, nm_pdf))
  }
}

sims <- sims %>% 
  bind_rows(.id = "country")

# Rename the age groups
age_nm <- levels(sims$age_cat)
age_nm <- str_extract_all(string = age_nm, pattern = "\\d+\\.\\d+|\\d+", simplify = T)

sims <- sims %>% 
  mutate(country = fct_recode(country, 
                              "Czechia" = "Czech",
                              "UK" = "United-Kingdom", 
                              "USA" = "United_States"))

saveRDS(object = sims %>% filter(time >= max(time) - 19), 
        file = paste0("_saved/", nm_dir, "/all-", nm_dir, ".rds"))

# Subset simulated data -------------------------------------------------------------
age_select <- c("[20,40)", "[40,60)", "[60,Inf]")

sims_cur <- sims %>% 
  filter(time >= max(time) - 19, 
         var_nm %in% c("seroPrev", "seroPPV", "seroInc", "trueInc"), 
         age_cat %in% age_select) %>% 
  mutate(prop = n / pop)

# Summary: median of all variables
sims_sumry <- sims_cur %>% 
  group_by(country, var_nm, age_cat) %>% 
  summarise(med_prop = median(prop), 
            med_n = median(n)) %>% 
  ungroup()

# Print range of median Sp estimates
tmp <- sims_cur %>% 
  filter(var_nm == "seroPrev") %>% 
  group_by(age_cat) %>% 
  summarise(
    q_50 = quantile(prop, probs = 0.5), 
    q_inf = quantile(prop, probs = 0.025), 
    q_sup = quantile(prop, probs = 0.975)) %>% 
  ungroup()

print("Seroprevalence, 95% prediction interval")
print(tmp %>% mutate(q_inf = round(100 * q_inf, 1), q_sup = round(100 * q_sup, 1), q_50 = round(100 * q_50, 1)))

# Print range of median PPV estimates
tmp <- sims_cur %>% 
  filter(var_nm == "seroPPV") %>% 
  group_by(age_cat) %>% 
  summarise(
    q_50 = quantile(n, probs = 0.5, na.rm = T), 
    q_inf = quantile(n, probs = 0.025, na.rm = T), 
    q_sup = quantile(n, probs = 0.975, na.rm = T)) %>% 
  ungroup()

print("PPV, 95% prediction interval")
print(tmp %>% mutate(q_inf = round(100 * q_inf, 0), q_sup = round(100 * q_sup, 0), q_50 = round(100 * q_50, 0)))

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

pl <- ggplot(data = tmp %>% filter(country %in% levels(tmp$country)), 
             mapping = aes(x = 1e2 * val, y = country, fill = age_cat, color = age_cat)) + 
  geom_density_ridges(scale = 0.95,
                      stat = "density_ridges",
                      quantile_lines = T,
                      quantiles = 2,
                      jittered_points = add_rugs,
                      point_shape = "|",
                      position = position_points_jitter(height = 0),
                      alpha = 0.5) +
  # geom_point(alpha = 0.5) + 
  #stat_summary(fun.data = "median_hilow") + 
  #geom_boxplot() + 
  facet_wrap(~ var_nm, scales = "free_x", ncol = 2) + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  scale_fill_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "Value (%)", y = "Country", fill = "Age group (yr)", color = "") +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        legend.position = "top", 
        strip.text = element_text(size = rel(1))
  ) + 
  guides(color = F)
print(pl)

ggsave(filename = sprintf("_figures/%s/main_Sp_PPV.pdf", nm_dir), 
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

ggsave(filename = sprintf("_figures/%s/sup_Seroprev_Seroinc.pdf", nm_dir), 
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
             mapping = aes(x = age_cat, y = 1e5 * inc_mean, group = interaction(country, .id))) + 
  geom_line(alpha = 0.5, color = "grey") + 
  #scale_color_brewer(palette = "Paired") + 
  facet_wrap(~ var_nm, nrow = 2, scales = "free") + 
  theme_classic() + theme(strip.background = element_blank(), 
                          strip.text = element_text(size = rel(1)),
                          legend.position = c(0.8, 0.8)) + 
  #scale_y_sqrt() + 
  labs(x = "Age group", y = "Incidence rate (per year per 100,000)", color = "")
print(pl)

ggsave(filename = sprintf("_figures/%s/sup_age_distribution.pdf", nm_dir), 
       plot = pl, width = 8, height = 8)

# Sup figure: Sp breakdown: Vp, Rp1, and Rp2 -------------------------------------
tmp <- sims %>% 
  filter(time >= max(time) - 19, 
         var_nm %in% c("Rp1", "Rp2", "Vp"), 
         age_cat %in% age_select) %>% 
  mutate(prop = n / pop, 
         var_nm = factor(var_nm, levels = c("Vp", "Rp2", "Rp1"))) %>% 
  group_by(country, age_cat, var_nm) %>% 
  summarise(prop = mean(prop)) %>% 
  ungroup()

pl <- ggplot(data = tmp %>% filter(country == "UK"), 
             mapping = aes(x = age_cat, y = prop, fill = var_nm)) + 
  geom_col(position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = c("#fbb4ae", "#fed9a6", "#b3cde3"), 
                    labels = c("Immune boost, vaccinated (Vp)", "Immune boost, recovered (Rp2)", "True infection (Rp1)")) + 
  theme(legend.position = "top") + 
  #scale_y_sqrt() + 
  labs(x = "Age group", y = "Proportion (relative to seroprevalence)", fill = "") 
print(pl)

ggsave(filename = sprintf("_figures/%s/sup_Sp_breakdown.pdf", nm_dir), 
       plot = pl, width = 8, height = 8)


#######################################################################################################
# End
#######################################################################################################