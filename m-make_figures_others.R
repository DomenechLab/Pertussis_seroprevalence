#######################################################################################################
# Make other figures
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
debug_bool <- F
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)


# Plot demographic data ---------------------------------------------------
ct_list <- read.table(file = "_data/list_countries.txt") %>% pluck(1)
ct_list2 <- read.table(file = "_data/list_countries2.txt") %>% pluck(1)

ct_list <- c(ct_list, ct_list2)
ct_list <- sort(unique(ct_list))
rm(ct_list2)

demog <- vector(mode = "list", length = length(ct_list))
names(demog) <- ct_list

for(ct in ct_list) {
  demog[[ct]] <- read_csv(file = sprintf("_data/_demog/_2010/%s_country_level_age_distribution_85.csv", ct), 
                          col_names = c("age", "pop"), 
                          col_types = "d") %>% 
    arrange(age) %>% 
    filter(age <= 79)
}

names(demog)[names(demog) %in% c("United_States", "United-Kingdom", "Czech")] <- c("USA", "UK", "Czechia")

demog <- demog %>% 
  bind_rows(.id = "country") %>% 
  arrange(country) %>% 
  group_by(country) %>% 
  mutate(pop_rel = pop / sum(pop)) %>% 
  ungroup()

pl <- ggplot(data = demog, mapping = aes(x = age, y = 1e2 * pop_rel)) + 
  #geom_line() + 
  geom_col() + 
  facet_wrap(~ country, ncol = 2, scales = "fixed", axis.labels = "margins", axes = "all") + 
  theme_classic() + 
  theme(strip.background = element_blank()) + 
  labs(x = "Age (years)", y = "Population size (relative to total population size), %")
print(pl)

ggsave(plot = pl, 
       filename = "_figures/main/fig_demog_structure.pdf", width = 8, height = 8)

# Plot DTP3 coverage data ------------------------------------------
dtp3 <- read_xlsx(path = "_data/_vaccine_coverage/Diphtheria tetanus toxoid and pertussis (DTP) vaccination coverage 2024-11-09 10-33 UTC.xlsx")

dtp3 <- dtp3 %>% 
  select(-c(GROUP, ANTIGEN_DESCRIPTION, ANTIGEN)) %>% 
  filter(!is.na(NAME))

pl <- ggplot(data = dtp3, mapping = aes(x = YEAR, y = COVERAGE, color = COVERAGE_CATEGORY)) + 
  geom_line() + 
  facet_wrap(~ NAME, scales = "fixed", ncol = 3) + 
  labs(x = "Year", y = "DTP3 coverage (%)")
print(pl)

# Plot DTP4 coverage data -------------------------------------------------
dtp4 <- read_xlsx(path = "_data/_vaccine_coverage/Diphtheria tetanus toxoid and pertussis booster vaccination coverage 2024-16-09 12-05 UTC.xlsx")

dtp4 <- dtp4 %>% 
  select(-c(GROUP, ANTIGEN_DESCRIPTION, ANTIGEN)) %>% 
  filter(!is.na(NAME))

pl <- ggplot(data = dtp4, mapping = aes(x = YEAR, y = COVERAGE, color = COVERAGE_CATEGORY)) + 
  geom_line() + 
  facet_wrap(~ NAME, scales = "fixed", ncol = 3) + 
  labs(x = "Year", y = "DTP4 coverage (%)")
print(pl)

# Merge datasets ----------------------------------------------------------

dtp_all <- dtp3 %>% 
  filter(COVERAGE_CATEGORY == "WUENIC") %>% 
  select(-matches("COVERAGE_|TARGET|DOSES")) %>% 
  mutate(vaccine = "DTP3") %>% 
  full_join(y = dtp4 %>% 
              filter(COVERAGE_CATEGORY == "ADMIN", NAME != "Belgium") %>% 
              select(-matches("COVERAGE_|TARGET|DOSES")) %>% 
              mutate(vaccine = "DTP4"))

dtp_sumry <- dtp_all %>% 
  group_by(CODE, NAME, vaccine) %>% 
  filter(YEAR >= 1990) %>% 
  summarise(mean_cov = mean(COVERAGE, na.rm = T)) %>% 
  ungroup()

# Plot
pl <- ggplot(data = dtp_all %>% filter(YEAR >= 1990), 
             mapping = aes(x = YEAR, y = COVERAGE, linetype = vaccine)) + 
  geom_line() + 
  facet_wrap(~ NAME, scales = "fixed", ncol = 3) + 
  geom_text(data = dtp_sumry %>% filter(vaccine == "DTP3"), 
            mapping = aes(x = 2000, y = 50, label = paste0("E(v3) = ", round(mean_cov, 1), "%"))) + 
  geom_text(data = dtp_sumry %>% filter(vaccine == "DTP4"), 
            mapping = aes(x = 2000, y = 25, label = paste0("E(v4) = ", round(mean_cov, 1), "%"))) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(x = "Year", y = "DTP coverage (%)", linetype = "")
print(pl)


ggsave(filename = "_figures/_others/DTP_coverage.pdf", plot = pl, width = 10, height = 8)


#######################################################################################################
# END
#######################################################################################################