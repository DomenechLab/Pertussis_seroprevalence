#######################################################################################################
# Make figure for all scenarios (different rho values) across countries
# Pick the same value of D (duration of immunity) to facilitate the comparison
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

# Figure: Sp across countries, for different rho values (point size: SD) -------------------------------------------------------------
appender <- function(string) TeX(paste0("$\\rho = $", string))

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
  labs(x = "Seroprevalence (%)", y = "Country", color = "Age group (yrs)", size = "SD (%)")
print(pl)

ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-1.pdf", width = 8, height = 8)

# Figure: Sp across countries, for different rho values (color: PPV) -------------------------------------------------------------
appender <- function(string) TeX(paste0("$\\rho = $", string))

pl <- ggplot(data = sims_95_PI_wide, 
                 mapping = aes(x = 100 * Sp, y = fct_rev(country), shape = age_cat, color = 100 * PPV)) + 
  geom_point(size = rel(2)) + 
  facet_rep_wrap(~as.character(rho_val), 
                 labeller = as_labeller(x = appender, default = label_parsed), 
                 scales = "fixed") + 
  scale_color_viridis(option = "magma") + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country", shape = "Age group (yrs)", color = "PPV (%)")
print(pl)
ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-2.pdf", width = 8, height = 8)

# Figure: Sp across countries, for different rho values (point size: PPV) -------------------------------------------------------------
appender <- function(string) TeX(paste0("$\\rho = $", string))

pl <- ggplot(data = sims_95_PI_wide, 
             mapping = aes(x = 100 * Sp, y = fct_rev(country), color = age_cat, size = 100 * PPV)) + 
  geom_point() + 
  facet_rep_wrap(~as.character(rho_val), 
                 labeller = as_labeller(x = appender, default = label_parsed), 
                 scales = "fixed") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1) + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country", color = "Age group (yrs)", size = "PPV (%)")
print(pl)
ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-3.pdf", width = 8, height = 8)

# Figure: Time plot of Sp (convergence check) -------------------------------------------------

ct <- "USA"
tmp <- sims_wide %>% 
  filter(country == ct) %>% 
  mutate(rho_val = as.character(rho_val), 
         age_cat = factor(age_cat), 
         rho_val = factor(rho_val))
levels(tmp$age_cat) <- paste0(levels(tmp$age_cat), " yo")
levels(tmp$rho_val) <- paste0("rho = ", levels(tmp$rho_val))

# Plot overall incidence
pl <- ggplot(data = tmp, 
             mapping = aes(x = time, y = 1e2 * Sp, group = .id)) + 
  geom_line(color = "grey") + 
  facet_rep_grid(rho_val ~ age_cat, scales = "free_y") + 
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank()) +
  labs(x = "Time (years)", y = "Seroprevalence (%)")
print(pl)

ggsave(plot = pl, filename = sprintf("_figures/main/fig_Sp_convergence_%s.pdf", ct), width = 10, height = 8)

# Figure: Breakdown of Sp, for different rho values -------------------------------------------------
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
  labs(x = "Age group (yo)", y = "Relative proportion", fill = "") 
print(pl)
ggsave(plot = pl, filename = sprintf("_figures/main/fig_Sp_breakdown_%s.pdf", ct), width = 8, height = 8)

# Figure: (S1+S2), V, and R as a function of age and rho ----------------------------------------------------------------
ct <- "USA"

tmp <- sims_wide %>% 
  filter(country == ct) %>% 
  mutate(S = S1 + S2, 
         V_R = V + R, 
         Vp_Rp2 = Rp2 + Vp, 
         lambda = trueInc / S) %>% 
  select(rho_val:pop, S, Rp1, V_R, Vp_Rp2, seroPrev, lambda) %>% 
  pivot_longer(cols = S:lambda, names_to = "var_nm", values_to = "prop") %>% 
  mutate(prop = if_else(var_nm %in% c('lambda'), prop, prop / pop), 
         var_nm = factor(var_nm)) %>% 
  group_by(rho_val, age_cat, var_nm) %>% 
  summarise(prop = median(prop)) %>% 
  ungroup()

tmp <- tmp %>% 
  mutate(var_nm = fct_relevel(var_nm, "lambda", "seroPrev", "S", "Rp1", "V_R", "Vp_Rp2")) %>% 
  mutate(var_nm = fct_recode(var_nm, 
                             "Force~of~Infection~(lambda)" = "lambda",
                             "Fraction~susceptible~to~infection~(S[1]+S[2])" = "S",
                             "Seroprevalence~from~true~infections~(R[P1])" = "Rp1",
                             "Fraction~immune~to~infection~(V+R)" = "V_R",
                             "Seroprevalence~from~immune~boosts~(V[P]+R[P2])" = "Vp_Rp2",
                             "Overall~seroprevalence" = "seroPrev"
  ))

pl <- ggplot(data = tmp, 
             mapping = aes(x = age_cat, y = 100 * prop, color = rho_val, group = rho_val)) + 
  geom_point() + 
  geom_line() + 
  facet_rep_wrap(~ var_nm, scales = "free_y", ncol = 2, dir = "h", labeller = "label_parsed") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1) + 
  theme_classic() + 
  theme(legend.position = "top", 
        #panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Age group (yrs)", y = "Value (%)", color = TeX("Boosting coefficient ($\\rho$)"))
print(pl)

ggsave(plot = pl, filename = sprintf("_figures/main/fig_FoI_S_Sp_%s.pdf", ct), width = 8, height = 8)

#######################################################################################################
# End
#######################################################################################################