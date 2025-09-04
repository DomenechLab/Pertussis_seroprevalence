#######################################################################################################
# Make figure for all scenarios (different rho values) across countries
# Pick the same value of D (duration of immunity) to facilitate the comparison
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
library(lemon)
library(egg) # For function tag_facet
library(latex2exp)
library(ggnewscale) # for function new_scale_fill
library(colorspace)
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

levels(sims$country) <- sort(levels(sims$country))

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
  facet_wrap(~as.character(rho_val), 
             labeller = as_labeller(x = appender, default = label_parsed), 
             axes = "all", axis.labels = "margins",
             scales = "fixed") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1) + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country", color = "Age (years)", size = "SD (%)")
print(pl)

ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-1.pdf", width = 8, height = 8)

# Figure: Sp across countries, for different rho values (color: PPV) -------------------------------------------------------------
appender <- function(string) TeX(paste0("$\\rho = $", string))

pl <- ggplot(data = sims_95_PI_wide, 
             mapping = aes(x = 100 * Sp, y = fct_rev(country), shape = age_cat, color = 100 * PPV)) + 
  geom_point(size = rel(2)) + 
  facet_wrap(~as.character(rho_val), 
             labeller = as_labeller(x = appender, default = label_parsed), 
             axes = "all", axis.labels = "margins",
             scales = "fixed") + 
  scale_color_viridis(option = "magma") + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country", shape = "Age (years)", color = "PPV (%)")
print(pl)
ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-2.pdf", width = 8, height = 8)

# Figure: Sp across countries, for different rho values (point size: PPV) -------------------------------------------------------------
appender <- function(string) TeX(paste0("$\\rho = $", string))

pl <- ggplot(data = sims_95_PI_wide, 
             mapping = aes(x = 100 * Sp, y = fct_rev(country), color = age_cat, size = 100 * PPV)) + 
  geom_point() + 
  facet_wrap(~as.character(rho_val), 
             labeller = as_labeller(x = appender, default = label_parsed), 
             scales = "fixed", axes = "all", axis.labels = "margins") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1) + 
  theme_classic() + 
  theme(legend.position = "top", 
        panel.grid.major.y = element_line(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country SCM", color = "Age (years)", size = "PPV (%)")
print(pl)
ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-3.pdf", width = 8, height = 8)


# Figure: Sp across countries for different rho (color: PPV, one color scale per age group) ---------------------------------------------
appender <- function(string) TeX(paste0("$\\rho = $", string))
ages_order <- c("20-39", "40-59", "60-79")
ages_order <- ages_order[c(1,2,3)]
ages_order_nm <- paste0(ages_order, " yo")

pl <- ggplot(data = sims_95_PI_wide %>% filter(age_cat == ages_order[1]), 
             mapping = aes(x = 100 * Sp, y = fct_rev(country), color = 100 * PPV)) + 
  geom_point(size = rel(2)) + 
  facet_wrap(~as.character(rho_val), 
             labeller = as_labeller(x = appender, default = label_parsed), 
             axes = "all", axis.labels = "margins",
             scales = "fixed") + 
  scale_color_continuous_sequential(palette = "Purple-Blue", 
                                    name = ages_order_nm[1], 
                                    begin = 0.2, 
                                    limits = range(100 * sims_95_PI_wide$PPV[sims_95_PI_wide$age_cat == ages_order[1]])) +
  new_scale_color() + 
  geom_point(data = sims_95_PI_wide %>% filter(age_cat == ages_order[2]), 
             mapping = aes(x = 100 * Sp, y = fct_rev(country), color = 100 * PPV), 
             size = rel(2)) + 
  scale_color_continuous_sequential(palette = "Red-Yellow", 
                                    name = ages_order_nm[2], 
                                    begin = 0.2, 
                                    limits = range(100 * sims_95_PI_wide$PPV[sims_95_PI_wide$age_cat == ages_order[2]])) + 
  new_scale_color() + 
  geom_point(data = sims_95_PI_wide %>% filter(age_cat == ages_order[3]), 
             mapping = aes(x = 100 * Sp, y = fct_rev(country), color = 100 * PPV), 
             size = rel(2)) + 
  scale_color_continuous_sequential(palette = "Green-Yellow", 
                                    name = ages_order_nm[3],
                                    begin = 0.2, 
                                    limits = range(100 * sims_95_PI_wide$PPV[sims_95_PI_wide$age_cat == ages_order[3]])) + 
  theme_classic() + 
  theme(legend.position = "top", 
        plot.title = element_text(size = 11, hjust = 0.5),  
        panel.grid.major.y = element_line(),
        #panel.grid.major.y = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_text(size = 11)) +
  labs(x = "Seroprevalence (%)", y = "Country", title = "PPV (%)")
print(pl)
ggsave(plot = pl, filename = "_figures/main/fig_Sp_predictions_all_scenarios-4.pdf", width = 8, height = 8)

# Summary statistics ------------------------------------------------------
tmp <- sims_95_PI_wide %>% 
  filter(age_cat != "40-59") %>% 
  arrange(rho_val, country, age_cat) %>% 
  group_by(rho_val, country) %>% 
  summarise(Sp_rel = Sp[age_cat == "60-79"] / Sp[age_cat == "20-39"]) %>% 
  ungroup()

library(lme4)
mod_boost <- lmer(formula = log(PPV) ~ 1 + factor(rho_val) + (1 | country), 
                data = filter(sims_95_PI_wide, age_cat == "20-39"))

mod_age <- lmer(formula = log(PPV) ~ 1 + age_cat + (1 | country), 
                  data = filter(sims_95_PI_wide, rho_val == 0.5))

# Figure: Time plot of Sp (convergence check) -------------------------------------------------

ct <- "USA"
tmp <- sims_wide %>% 
  filter(country == ct) %>% 
  mutate(rho_val = as.character(rho_val), 
         age_cat = factor(age_cat), 
         rho_val = factor(rho_val))
levels(tmp$age_cat) <- paste0(levels(tmp$age_cat), " yo")
levels(tmp$rho_val) <- paste0("rho = ", levels(tmp$rho_val))

tmp <- tmp %>% 
  mutate(seroInc = seroInc / pop, 
         trueInc = trueInc / pop) %>% 
  pivot_longer(cols = -c(rho_val:pop))

# Plot overall incidence
pl <- ggplot(data = tmp %>% filter(name %in% c("seroInc", "trueInc")), 
             mapping = aes(x = time, y = 1e5 * value, 
                           group = interaction(.id, name), color = name)) + 
  geom_line() + 
  facet_rep_grid(rho_val ~ age_cat, scales = "free_y") + 
  theme_classic() + 
  theme(legend.position = "top", 
        strip.background = element_blank()) +
  scale_y_log10() + 
  scale_color_manual(labels = c("Sero-incidence", "True incidence of infections"), 
                     values = c("seroInc" = "#fb8072", "trueInc" = "#80b1d3")) + 
  labs(x = "Time (years)", y = "Incidence rate (per 100,000)", color = "")
print(pl)

ggsave(plot = pl, filename = sprintf("_figures/main/fig_Sp_convergence_%s.pdf", ct), 
       width = 10, height = 8)

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
  labs(x = "Age (years)", y = "Relative proportion", fill = "") 
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
  mutate(var_nm2 = fct_recode(var_nm, 
                              "Force~of~Infection~(lambda)" = "lambda",
                              "Fraction~susceptible~to~infection~(S[1]+S[2])" = "S",
                              "Seroprevalence~from~true~infections~(R[P1])" = "Rp1",
                              "Fraction~immune~to~infection~(V+R)" = "V_R",
                              "Seroprevalence~from~immune~boosts~(V[P]+R[P2])" = "Vp_Rp2",
                              "Overall~seroprevalence" = "seroPrev"), 
         var_nm3 = fct_recode(var_nm, 
                              "lambda" = "lambda",
                              "S[1] + S[2]" = "S",
                              "R[P1]" = "Rp1",
                              "R + V" = "V_R",
                              "R[P2] + V[P]" = "Vp_Rp2",
                              "R[P1] +  R[P2] + V[P]" = "seroPrev"))

pl <- ggplot(data = tmp, 
             mapping = aes(x = age_cat, y = 100 * prop, color = rho_val, group = rho_val)) + 
  geom_point() + 
  geom_line() + 
  # geom_text(data = data.frame(text = LETTERS[1:6], var_nm3 = levels(tmp$var_nm3)), 
  #           mapping = aes(x = -Inf, y = Inf, label = text), 
  #           hjust = -0.2, vjust = -0.2, 
  #           inherit.aes = F) + 
  facet_rep_wrap(~ var_nm3, scales = "free_y", ncol = 2, 
                 dir = "h", labeller = "label_parsed", strip.position = "left") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1)

pl <- tag_facet(p = pl, open = "", close = "", tag_pool = LETTERS, fontface = 2)

pl <- pl + 
  theme_classic() + 
  #theme_bw() + 
  theme(legend.position = "top", 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(), 
        strip.placement = "outside",
        #strip.background = element_rect(fill = "#f0f0f0"),
        strip.text = element_text(size = 11)) +
  labs(x = "Age (years)", y = "", color = TeX("Boosting coefficient ($\\rho$)"))
print(pl)

ggsave(plot = pl, filename = sprintf("_figures/main/fig_FoI_S_Sp_%s.pdf", ct), 
       width = 8, height = 8)

# Same plot, other display ------------------------------------------------
# df_nm <- data.frame(var_nm = levels(tmp$var_nm),
#                     var_nm2 = c("Force of infection, $\\lambda$ (% per year)",
#                                 "Fraction susceptible to infection, $S_1 + S_2$ (%)",
#                                 "Seroprevalence from true infections, $R_{P,1}$ (%)",
#                                 "Fraction immune to infection, V+R (%)",
#                                 "Seroprevalence from immune boosts, $V_P + R_{P,2}$ (%)",
#                                 "Overall seroprevalence"), 
#                     var_nm3 = c("$\\lambda$ (% per year)",
#                                 "$S_1 + S_2$ (%)",
#                                 "$R_{P,1}$ (%)",
#                                 "$V+R$ (%)",
#                                 "$R_{P,2} + V_P$ (%)",
#                                 "$R_{P,1} + R_{P,2} + V_P$ (%)"))
# 
# pl <- vector(mode = "list", length = nlevels(tmp$var_nm))
# for(i in seq_along(levels(tmp$var_nm))) {
#   
#   var_cur <- levels(tmp$var_nm)[i]
#   pl[[i]] <- ggplot(data = tmp %>% filter(var_nm == var_cur),
#                     mapping = aes(x = age_cat, y = 100 * prop, color = rho_val, group = rho_val)) +
#     geom_point() +
#     geom_line() +
#     scale_color_viridis(option = "viridis", discrete = T, end = 1, direction = -1) +
#     theme_classic() +
#     theme(legend.position = "top",
#           #panel.grid.major.y = element_line(),
#           strip.background = element_blank(),
#           strip.text = element_text(size = 11)) +
#     labs(x = "Age group (years)", y = TeX(df_nm$var_nm3[i]), 
#          #color = TeX("Boosting coefficient ($\\rho$)")
#          color = "Boosting coefficient"
#     )
#   print(pl[[i]])
# }
# 
# # Arrange plots
# 
# pl_all <- (pl[[1]] / pl[[3]] / pl[[5]] + plot_layout(axes = "collect_x")) | 
#   (pl[[2]] / pl[[4]] / pl[[6]] + plot_layout(axes = "collect_x"))
# pl_all <- pl_all +  plot_layout(guides = "collect", axes = "collect") & theme(legend.position='top') 
# pl_all <- pl_all + plot_annotation(tag_levels = "A") 
# print(pl_all)
# 
# ggsave(plot = pl_all, filename = sprintf("_figures/main/fig_FoI_S_Sp_%s.pdf", ct), 
#        width = 8, height = 8)

#######################################################################################################
# End
#######################################################################################################