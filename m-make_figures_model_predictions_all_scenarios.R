#######################################################################################################
# Make figure for all scenarios (different rho values) across countries
# Pick the best estimate of 1 / alpha (with lowest MWSE) for every value of rho
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
library(lemon)
library(latex2exp)
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)

# Load results ------------------------------------------------------------
rho_choose <- c(0.5, 1, 2, 5)
D_choose <- c(40, 40, 40, 40)

sims <- vector(mode = "list", length = length(rho_choose))
names(sims) <- as.character(rho_choose)

for(i in seq_along(rho_choose)) {
  dir_select <- sprintf("_saved/prediction_12_countries/DV_%d-rhoV_%.2f-DR_%d-rhoR_%.2f/", 
                        D_choose[i], rho_choose[i], D_choose[i], rho_choose[i])
  file_select <- list.files(path = dir_select, pattern = "all-")
  print(file_select)
  
  sims[[i]] <- readRDS(paste0(dir_select, file_select))
}

sims <- bind_rows(sims, .id = "rho_val")
age_cats <- levels(sims$age_cat)
age_cat_choose <- age_cats[(length(age_cats) - 2):length(age_cats)]
print(age_cat_choose)

# Recast in wide format and calculate PPV ---------------------------------
sims_wide <- sims %>% 
  filter(age_cat %in% age_cat_choose) %>% 
  select(-var_type) %>% 
  pivot_wider(names_from = "var_nm", values_from = "n") %>% 
  mutate(Sp = seroPrev / pop, 
         PPV = Rp1 / seroPrev, 
         age_cat = fct_recode(age_cat, 
                              "20-39" = "[20,40)", 
                              "40-59" = "[40,60)", 
                              "60-79" = "[60,Inf]"))

# Calculate 95% prediction intervals for Sp and PPV
sims_95_PI <- sims_wide %>% 
  select(rho_val:pop, Sp, PPV) %>% 
  pivot_longer(cols = c("Sp", "PPV"), names_to = "var_nm", values_to = "prop") %>% 
  group_by(rho_val, country, age_cat, var_nm) %>% 
  summarise(q_inf = quantile(prop, probs = 0.25, na.rm = T), 
            q_med = quantile(prop, probs = 0.5, na.rm = T), 
            q_sd = sd(prop, na.rm = T),
            q_sup = quantile(prop, probs = 0.75, na.rm = T)) %>% 
  ungroup()

sims_95_PI_wide <- sims_95_PI %>% 
  select(-c(q_inf, q_sup, q_sd)) %>% 
  pivot_wider(names_from = "var_nm", values_from = "q_med")

# Figure: Sp across countries, for different rho values -------------------------------------------------------------
appender <- function(string) TeX(paste0("$\\rho_R = \\rho_V = $", string))

pl <- ggplot(data = sims_95_PI %>% filter(var_nm == "Sp"), 
             mapping = aes(x = 100 * q_med, y = fct_rev(country), color = age_cat, size = 100 * q_sd)) + 
  geom_point() + 
  #geom_linerange(mapping = aes(xmin = 100 * q_inf, xmax = 100 * q_sup)) + 
  facet_rep_wrap(~as.character(rho_val), 
                 labeller = as_labeller(x = appender, default = label_parsed), 
                 scales = "fixed") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1) + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country", color = "Age group (yr)", size = "SD (%)")
print(pl)

# Figure: Breakdown of Sp, for different rho values -------------------------------------------------
ct <- "USA"

# Plot overall incidence
pl <- ggplot(data = sims_wide %>% filter(country == ct), 
             mapping = aes(x = time, y = 1e5 * trueInc / pop, color = .id)) + 
  geom_line() + 
  facet_grid(rho_val ~ age_cat, scales = "free")
print(pl)

tmp <- sims_wide %>% 
  filter(country == ct) %>% 
  select(rho_val:pop, Vp, Rp1, Rp2) %>% 
  pivot_longer(cols = Vp:Rp2, names_to = "var_nm", values_to = "prop") %>% 
  mutate(prop = prop / pop, 
         var_nm = factor(var_nm, levels = c("Vp", "Rp2", "Rp1"))) %>% 
  group_by(rho_val, age_cat, var_nm) %>% 
  summarise(prop = median(prop)) %>% 
  ungroup()

pl <- ggplot(data = tmp, 
             mapping = aes(x = age_cat, y = prop, fill = var_nm)) + 
  geom_col(position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = c("#fbb4ae", "#fed9a6", "#b3cde3"), 
                    labels = c("Immune boost, vaccinated (Vp)", "Immune boost, recovered (Rp2)", "True infection (Rp1)")) + 
  facet_rep_wrap(~as.character(rho_val), 
                 labeller = as_labeller(x = appender, default = label_parsed), 
                 scales = "fixed") + 
  theme(legend.position = "top", 
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Age group (yr)", y = "Relative proportion", fill = "") 
print(pl)

#######################################################################################################
# End
#######################################################################################################