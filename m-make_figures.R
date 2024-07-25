#######################################################################################################
# Make figures
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
add_rugs <- F # Should rugs be added for every data point

# Load simulations -------------------------------------------------------------
l_countries <- list.dirs(path = "_saved", full.names = F, recursive = F)
sims <- vector(mode = "list", length = length(l_countries))
names(sims) <- l_countries

for(i in seq_along(l_countries)) {
  sims[[i]] <- readRDS(file = sprintf("_saved/%s/alphaV_0.02-rhoV_0.5-vacCov_0.95-2doses.rds", l_countries[i])) %>% 
    pluck("merged_ages")
}

sims <- sims %>% 
  bind_rows(.id = "country") %>% 
  mutate(country = fct_recode(country, 
                              "S. Africa" = "South_Africa", 
                              "UK" = "United-Kingdom", 
                              "USA" = "United_States"))

# Make figure -------------------------------------------------------------
sims_cur <- sims %>% 
  filter(time >= max(time) - 19, 
         var_nm %in% c("seroPrev", "seroPPV", "seroInc", "trueInc"), 
         age_cat %in% c("[20,40)", "[40,60)", "[60,Inf]")) %>% 
  mutate(prop = n / pop)

# Summary: median of all variables
sims_sumry <- sims_cur %>% 
  group_by(country, var_nm, age_cat) %>% 
  summarise(med_prop = median(prop), 
            med_n = median(n)) %>% 
  ungroup()

# Order countries by increasing seroprevalence
country_order <- sims_sumry %>% 
  filter(var_nm == "seroPrev", age_cat == "[20,40)") %>% 
  arrange(med_prop) %>% 
  pluck("country")

sims_cur <- sims_cur %>% 
  mutate(country = factor(x = country, levels = country_order))


# Plot incidence vs. seroprevalence ---------------------------------------
sims_cur_long <- sims_cur %>% 
  select(-c(var_type, prop)) %>% 
  pivot_wider(names_from = "var_nm", values_from = "n")

pl1 <- ggplot(data = sims_cur_long, 
              mapping = aes(x = seroPrev / pop, y = seroInc / pop, color = country)) + 
  geom_point(alpha = 1) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  facet_wrap(~ age_cat, ncol = 2, scales = "fixed") + 
  scale_color_viridis(option = "magma", discrete = T) + 
  labs(x = "Seroprevalence", y = "Sero-incidence (per year)")
print(pl1)

pl2 <- ggplot(data = sims_cur_long, 
              mapping = aes(x = seroPPV, y = trueInc / seroInc, color = country)) + 
  geom_point(alpha = 1) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  facet_wrap(~ age_cat, ncol = 2, scales = "fixed") + 
  scale_color_viridis(option = "magma", discrete = T) + 
  labs(x = "PPV", y = "True incidence / sero-incidence")
print(pl2)

# Left figure -------------------------------------------------------------

tmp <- sims_cur %>% 
  filter(var_nm %in% c("seroPrev", "seroPPV")) %>% 
  mutate(val = ifelse(var_nm == "seroPrev", prop, n), 
         var_nm = factor(var_nm, levels = c("seroPrev", "seroPPV"), labels = c("Seroprevalence", "PPV of serology")))

pl_left <- ggplot(data = tmp, 
                  mapping = aes(x = 1e2 * val, y = country, fill = age_cat, color = age_cat)) + 
  geom_density_ridges(scale = 0.9,
                      stat = "density_ridges", 
                      quantile_lines = T, 
                      quantiles = 2, 
                      jittered_points = add_rugs,
                      point_shape = "|", 
                      position = position_points_jitter(height = 0),
                      alpha = 0.5) + 
  facet_wrap(~ var_nm, scales = "free_x", ncol = 1) + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  scale_fill_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "Value (%)", y = "Country", fill = "Age group", color = "") +
  theme_classic() + 
  theme(strip.background = element_blank(), legend.position = "top") + 
  guides(color = F)
print(pl_left)

# Seroprevalence 
# pl_prev <- ggplot(data = sims_cur %>% filter(var_nm == "seroPrev"), 
#                   mapping = aes(x = 1e2 * prop, y = country, fill = age_cat, color = age_cat)) + 
#   geom_density_ridges(scale = 0.9,
#                       stat = "density_ridges", 
#                       quantile_lines = T, 
#                       quantiles = 2, 
#                       jittered_points = add_rugs,
#                       point_shape = "|", 
#                       position = position_points_jitter(height = 0),
#                       alpha = 0.5) + 
#   scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
#   scale_fill_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
#   labs(x = "Seroprevalence (%)", y = "Country") +
#   theme_classic()
# print(pl_prev)
# 
# # PPV
# pl_ppv <- ggplot(data = sims_cur %>% filter(var_nm == "seroPPV"), 
#                  mapping = aes(x = 1e2 * n, y = country)) + 
#   geom_density_ridges(scale = 0.9,
#                       stat = "density_ridges", 
#                       quantile_lines = T, 
#                       quantiles = 2, 
#                       jittered_points = add_rugs,
#                       point_shape = "|", 
#                       position = position_points_jitter(height = 0),
#                       alpha = 0.5) + 
#   labs(x = "Positive predictive value (%)", y = "Country") +
#   theme_classic()
# print(pl_ppv)


# Right plot: Incidence outcomes --------------------------------------------------------------
tmp <- sims_cur %>% 
  filter(var_nm %in% c("seroInc", "trueInc")) %>% 
  mutate(var_nm = factor(var_nm, 
                         levels = c("trueInc", "seroInc") %>% rev(), 
                         labels = c("True incidence", "Sero-incidence") %>% rev()), 
         country_var = interaction(country, var_nm, sep = "_", lex.order = T))

y_breaks <- levels(tmp$country_var) %>% str_split_i(pattern = "_", i = 1)
y_breaks[c(T, F)] <- ""

pl_right <- ggplot(data = tmp, 
                   mapping = aes(x = 1e5 * prop, y = country_var, color = age_cat, fill = age_cat, linetype = var_nm)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      jittered_points = add_rugs,
                      point_shape = "|", 
                      position = position_points_jitter(height = 0),
                      stat = "density_ridges", 
                      scale = 0.9,
                      alpha = 0.5) + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  scale_fill_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "Incidence rate (per year per 100,000)", 
       y = "", 
       linetype = "") +
  guides(color = F, fill = F) + 
  scale_x_log10() +
  scale_y_discrete(labels = y_breaks) + 
  theme_classic() + 
  theme(legend.position = "top")
print(pl_right)

pl_right_bis <- ggplot(data = tmp, 
                       mapping = aes(x = 1e5 * prop, y = country, color = age_cat, fill = age_cat)) + 
  geom_density_ridges(quantile_lines = T, 
                      quantiles = 2, 
                      jittered_points = add_rugs,
                      point_shape = "|", 
                      position = position_points_jitter(height = 0),
                      stat = "density_ridges", 
                      scale = 0.9,
                      alpha = 0.5) + 
  facet_wrap(~ var_nm, ncol = 1, scales = "free_x") + 
  scale_color_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  scale_fill_viridis(option = "viridis", discrete = T, end = 0.5, direction = -1) + 
  labs(x = "Incidence rate (per year per 100,000)", 
       y = "", 
       linetype = "") +
  guides(color = F, fill = F) + 
  #scale_x_log10() +
  theme_classic() + 
  theme(legend.position = "top", strip.background = element_blank())
print(pl_right_bis)

# Piece the plots together ------------------------------------------------
pl_all <- (pl_left | pl_right_bis) + 
  plot_layout(axes = "collect", ncol = 2, guides = "auto") & 
  theme(legend.position = "top")

print(pl_all)
ggsave(filename = "_figures/test3.pdf", plot = pl_all, width = 10, height = 8)

#######################################################################################################
# End
#######################################################################################################