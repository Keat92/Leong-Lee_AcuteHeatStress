# =============================================================================
# Part 3 — Colony ED50 vs field bleaching (2024 heatwave)
# =============================================================================
# Links each colony's assay ED50 (Part 2) to the field monitoring of the same
# tagged colonies (May–November 2024) and asks whether more heat-tolerant
# colonies (higher ED50) bleached less in the field:
#   * all ED50 comparisons use the July–August surveys (peak bleaching):
#     Spearman correlation of peak bleaching ~ ED50, GLM / GLMM
#     of bleaching ~ ED50 * species, and bleaching by within-species ED50
#     group (Fisher + Cochran-Armitage),
#   * descriptive: survey coverage and seasonal bleaching index (BI) per
#     species and site, May–November.
#
# Moved from Part6_bleaching.R. The unused bleachingindex.csv import was
# dropped (the seasonal BI is recomputed below from the raw field sheet).
#
# Inputs:
#   data/processed/colony_ed50_ranked.csv  colony ED50 + within-species rank
#                                          (written by part2_drc_ed50analysis.R)
#   data/raw/field_monitoring.csv          field surveys of tagged colonies
#                                          (sheet "Fieldmonitoringdatasheet" of
#                                          Fieldmonitoringdatasheet.xlsx)
# Outputs:
#   data/processed/field_ed50_compare.csv  field surveys joined to colony ED50
#   output/tables/   Spearman, Fisher, survey coverage, bleaching-% and BI tables
#   output/figures/  ED50 vs bleaching, quartile and BI figures
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)  # combine panels
library(rstatix)    # cor_test(), shapiro_test(), wilcox_test()
library(lme4)       # glmer()
library(broom)      # tidy()
library(car)        # Anova()
library(here)
# MASS::glmmPQL() is called with `::` so MASS::select() never masks dplyr.

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# Colony ED50 from Part 2. Quartiles are assigned within species, so every
# species is split by its own ED50 distribution.
# Q1 = highest ED50 (most heat tolerant), Q4 = lowest ED50 (least tolerant).
# cut() bins from low to high, so labels run Q4 -> Q1; levels are then
# reordered Q1 -> Q4. `include.lowest = TRUE` keeps the minimum colony in Q4.
ed50 <- read_csv(here("data", "processed", "colony_ed50_ranked.csv"),
                 show_col_types = FALSE) %>%
  group_by(Species) %>%
  mutate(Quartile = cut(ED50_mean,
                        breaks = quantile(ED50_mean, probs = seq(0, 1, 0.25),
                                          na.rm = TRUE),
                        labels = c("Q4", "Q3", "Q2", "Q1"),
                        include.lowest = TRUE),
         Quartile = factor(Quartile, levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
  ungroup()
str(ed50)
table(ed50$Species, ed50$Quartile)

# Field surveys. Bleach% is read as text because "-" marks "not scored".
month_order <- c("May", "June", "July", "August", "September", "November")

field_raw <- read_csv(here("data", "raw", "field_monitoring.csv"),
                      col_types = cols(`Bleach%` = col_character()))
str(field_raw)

# -----------------------------------------------------------------------------
# 2. Field monitoring of tagged colonies
# -----------------------------------------------------------------------------
# keep survey columns; "-" in Bleach% means not scored -> NA.
# Bleach% and Mortality% are converted to proportions (0-1) here, straight from
# the raw percentages, so re-running later lines can never divide by 100 twice.
field <- field_raw %>%
  select(month = Month, location = Location, species = Species,
         colony = Colony, colour = Colour,
         bleach = `Bleach%`, mortality = `Mortality%`) %>%
  mutate(month = factor(month, levels = month_order),
         bleach = as.numeric(na_if(bleach, "-")) / 100,
         mortality = mortality / 100)

# standardise colony IDs: drop spaces/"?" and zero-pad number (GOS1 -> GOS01)
fix_colony <- function(x) {
  x <- gsub("[^A-Za-z0-9]", "", x)
  paste0(sub("[0-9]+$", "", x),
         sprintf("%02d", as.integer(sub("^[A-Za-z]+", "", x))))
}

# add colony ED50, quartile and rank to field data
field <- field %>%
  mutate(colony = fix_colony(colony)) %>%
  left_join(ed50 %>%
              mutate(colony = fix_colony(Colony)) %>%
              select(colony, ed50_mean = ED50_mean,
                     quartile = Quartile, rank = Rank),
            by = "colony")
str(field)

# field colonies without an ED50 value
field %>% filter(is.na(ed50_mean)) %>% distinct(location, colony)

# keep only colonies with an ED50 value
field <- field %>% filter(!is.na(ed50_mean))

# bleach is a proportion (0-1); lower = less bleaching.
# (mortality is kept in field_ed50_compare.csv but not analysed here)
summary(field$bleach)

# bleaching status: 0 = bleaching (bleach > 0), 1 = healthy (bleach = 0)
field <- field %>%
  mutate(bleaching_occurrence = ifelse(bleach > 0, 0, 1))
table(field$month, field$bleaching_occurrence, useNA = "ifany")

write_csv(field, here("data", "processed", "field_ed50_compare.csv"))

# -----------------------------------------------------------------------------
# 3. July–August subset (peak bleaching period)
# -----------------------------------------------------------------------------
# All ED50–field comparisons below use the July and August surveys only.
peak_sub <- field %>%
  filter(month %in% c("July", "August")) %>%
  mutate(month = droplevels(month))
str(peak_sub)
count(peak_sub, month, species)

field_peak_ja <- peak_sub %>%
  filter(!is.na(bleach)) %>%
  group_by(species, location, colony, ed50_mean, quartile) %>%
  summarise(bleach_max = max(bleach), n_surveys = n(), .groups = "drop") %>%
  mutate(bleaching_occurrence = ifelse(bleach_max > 0, 0, 1))  # 0 = bleaching, 1 = healthy
count(field_peak_ja, species)
count(field_peak_ja, n_surveys)

# -----------------------------------------------------------------------------
# 4. Spearman correlation: colony ED50 vs peak bleaching (Jul–Aug)
# -----------------------------------------------------------------------------
# Shapiro-Wilk normality test on peak bleach proportion (all species pooled,
# then per species). Mostly zeros (not normal), so use Spearman.
field_peak_ja %>%
  shapiro_test(bleach_max)

field_peak_ja %>%
  group_by(species) %>%
  shapiro_test(bleach_max)

# within each species
cor_ja <- field_peak_ja %>%
  group_by(species) %>%
  cor_test(ed50_mean, bleach_max, method = "spearman")
cor_ja

# all species pooled
cor_ja_all <- field_peak_ja %>%
  cor_test(ed50_mean, bleach_max, method = "spearman") %>%
  mutate(species = "All species")
cor_ja_all

### Scatterplot: colony ED50 vs peak bleaching ###
species_cols <- c(GOS = "#2a78d6", PAC = "#eb6834", POC = "#1baf7a", POR = "#eda100")

# legend labels carry each species' Spearman result
spearman_labs <- with(cor_ja, setNames(
  sprintf("%s: rho = %.2f, p = %.3f", species, cor, p), species))

overall_lab <- with(cor_ja_all, sprintf(
  "All species: rho = %.2f, %s", cor,
  ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))))

fig_cor_ja <- ggplot(field_peak_ja, aes(ed50_mean, bleach_max, colour = species)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey20", linewidth = 0.8) +
  annotate("text", x = Inf, y = 1, hjust = 1, vjust = 1, size = 3.5,
           colour = "grey20", label = overall_lab) +
  scale_colour_manual(values = species_cols, labels = spearman_labs, name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Colony mean ED50 (°C)", y = "Peak bleaching proportion (Jul-Aug)") +
  guides(colour = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.justification = "left")
fig_cor_ja
ggsave(here("output", "figures", "fig_ed50_spearman_JulAug.png"), fig_cor_ja,
       width = 6, height = 5.5, dpi = 300)

### Table: Spearman correlation of colony ED50 vs peak bleaching ###
n_bleached <- bind_rows(field_peak_ja,
                        mutate(field_peak_ja, species = "All species")) %>%
  group_by(species) %>%
  summarise(n = n(), n_bleached = sum(bleach_max > 0), .groups = "drop")

cor_table <- bind_rows(cor_ja, cor_ja_all) %>%
  left_join(n_bleached, by = "species") %>%
  transmute(Species = species, n, n_bleached,
            rho = round(cor, 2), p = signif(p, 2))
cor_table
write_csv(cor_table, here("output", "tables", "spearman_ed50_bleach_JulAug.csv"))

# -----------------------------------------------------------------------------
# 5. Models: bleaching ~ ED50 * species (July + August surveys)
# -----------------------------------------------------------------------------
### Bleaching occurrence: 0 = bleaching, 1 = healthy; binomial GLM ###
m_occ_ja <- glm(bleaching_occurrence ~ ed50_mean * species,
                family = binomial, data = peak_sub)
summary(m_occ_ja)
anova(m_occ_ja, test = "Chisq")   # sequential (Type I)
drop1(m_occ_ja, test = "Chisq")   # interaction, adjusted
drop1(update(m_occ_ja, . ~ ed50_mean + species), test = "Chisq")  # main effects, adjusted

### Bleaching proportion ###
# bleach = proportion of colony area (0-1), not counts -> quasibinomial

## GLM (repeat surveys of the same colony treated as independent)
m_bleach_ja <- glm(bleach ~ ed50_mean * species,
                   family = quasibinomial, data = peak_sub)
summary(m_bleach_ja)
anova(m_bleach_ja, test = "F")    # sequential (Type I)
drop1(m_bleach_ja, test = "F")    # interaction, adjusted
drop1(update(m_bleach_ja, . ~ ed50_mean + species), test = "F")  # main effects, adjusted

## GLMM with colony as random effect (colonies surveyed in both months)
# glmer() has no quasi families -> glmmPQL (penalised quasi-likelihood)
m_bleach_ja_mixed <- MASS::glmmPQL(bleach ~ ed50_mean + species,
                                   random = ~ 1 | colony,
                                   family = quasibinomial, data = peak_sub,
                                   verbose = FALSE)
summary(m_bleach_ja_mixed)
car::Anova(m_bleach_ja_mixed, type = 2)   # Wald chi-square (anova() not available)

# -----------------------------------------------------------------------------
# 6. Peak bleaching (Jul–Aug) by ED50 group: Q1, Q2-Q3, Q4
# -----------------------------------------------------------------------------
# Q2 and Q3 pooled; Q1 = highest ED50 (high tolerance),
# Q4 = lowest ED50 (low tolerance)
tolerance_group <- function(q) {
  factor(case_when(q == "Q1" ~ "Q1",
                   q %in% c("Q2", "Q3") ~ "Q2-Q3",
                   q == "Q4" ~ "Q4"),
         levels = c("Q1", "Q2-Q3", "Q4"))
}

### ED50 range per group: Q1, Q2-Q3, Q4 (all colonies in the ED50 assay) ###
species_labels <- c(GOS = "G. pectinata", POC = "P. acuta",
                    POR = "P. lutea", PAC = "P. speciosa")
group3_cols <- c(Q1 = "#B2182B", "Q2-Q3" = "#9E9E9E", Q4 = "#2166AC")  # red = high ED50, as in Part5

group3_range <- ed50 %>%
  mutate(ed50_group = tolerance_group(Quartile)) %>%
  group_by(Species, ed50_group) %>%
  summarise(n_colony    = n(),
            ED50_min    = min(ED50_mean),
            ED50_max    = max(ED50_mean),
            ED50_median = median(ED50_mean),
            ED50_span   = sprintf("%.2f-%.2f", ED50_min, ED50_max),
            .groups = "drop")
group3_range
write_csv(group3_range, here("output", "tables", "colony_quartile3_range.csv"))

fig_group3_range <- group3_range %>%
  mutate(SpeciesLabel = factor(species_labels[Species],
                               levels = species_labels[sort(names(species_labels))])) %>%
  ggplot(aes(x = ed50_group, colour = ed50_group)) +
  geom_linerange(aes(ymin = ED50_min, ymax = ED50_max),
                 linewidth = 6, alpha = 0.6) +
  geom_point(aes(y = ED50_median), size = 2, colour = "black") +
  facet_wrap(~ SpeciesLabel, nrow = 1) +   # shared y-scale: aligns with panel a
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  scale_colour_manual(values = group3_cols) +
  labs(x = expression(ED[50]~quartile), y = expression(ED[50]~range~"(°C)")) +
  theme_classic() +
  theme(strip.text = element_text(face = "italic", size = 11),
        panel.spacing.x = unit(1, "lines"),
        legend.position = "none")
fig_group3_range
ggsave(here("output", "figures", "fig_quartile3_range.png"), fig_group3_range, width = 9, height = 3.5, dpi = 300)


### Stacked bar: bleached vs unbleached by ED50 group (July-August peak) ###
peak_ja_bleach <- field_peak_ja %>%
  mutate(ed50_group = tolerance_group(quartile)) %>%
  group_by(species, ed50_group) %>%
  summarise(n = n(), n_bleached = sum(bleach_max > 0), .groups = "drop")
peak_ja_bleach

peak_ja_stack <- peak_ja_bleach %>%
  mutate(Bleached = n_bleached / n, Unbleached = 1 - Bleached) %>%
  pivot_longer(c(Bleached, Unbleached),
               names_to = "status", values_to = "prop") %>%
  mutate(status = factor(status, levels = c("Unbleached", "Bleached")))

fig_peak_ja_stack <- ggplot(peak_ja_stack, aes(ed50_group, prop, fill = status)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.5) +
  geom_text(data = peak_ja_bleach, aes(ed50_group, 1, label = paste0("n = ", n)),
            inherit.aes = FALSE, vjust = -0.5, size = 3, colour = "grey30") +
  facet_wrap(~ species, nrow = 1, labeller = as_labeller(species_labels)) +
  scale_fill_manual(values = c(Unbleached = "grey80", Bleached = "grey25"),
                    breaks = c("Bleached", "Unbleached"), name = NULL) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression(ED[50]~quartile~"(Q1 = high tolerance, Q4 = low tolerance)"),
       y = "Proportion of colonies") +
  theme_classic() +
  theme(legend.position = "top",
        strip.text = element_text(face = "italic", size = 11),
        panel.spacing.x = unit(1.5, "lines"))
fig_peak_ja_stack
ggsave(here("output", "figures", "fig_quartile3_bleached_peak_JulAug.png"), fig_peak_ja_stack,
       width = 9, height = 4.5, dpi = 300)
ggsave(here("output", "figures", "fig_quartile3_bleached_peak_JulAug.pdf"), fig_peak_ja_stack,
       width = 9, height = 4.5)


### Fisher's exact + Cochran-Armitage trend: ED50 group vs bleached (July-August peak) ###
fisher_peak_ja <- field_peak_ja %>%
  mutate(ed50_group = tolerance_group(quartile),
         bleached = factor(bleach_max > 0, levels = c(TRUE, FALSE),
                           labels = c("Bleached", "Unbleached"))) %>%
  group_by(species) %>%
  group_modify(~ {
    tab <- table(.x$ed50_group, .x$bleached)              # 3 x 2 table per species
    tibble(n = sum(tab), n_bleached = sum(tab[, "Bleached"]),
           bleached_Q1   = paste0(tab["Q1", "Bleached"], "/", sum(tab["Q1", ])),
           bleached_Q2Q3 = paste0(tab["Q2-Q3", "Bleached"], "/", sum(tab["Q2-Q3", ])),
           bleached_Q4   = paste0(tab["Q4", "Bleached"], "/", sum(tab["Q4", ])),
           p       = fisher.test(tab)$p.value,
           p_trend = prop.trend.test(tab[, "Bleached"], rowSums(tab), score = 1:3)$p.value)
  }) %>%
  ungroup() %>%
  mutate(p_holm       = p.adjust(p, method = "holm"),
         p_trend_holm = p.adjust(p_trend, method = "holm"))
fisher_peak_ja
write_csv(fisher_peak_ja, here("output", "tables", "fisher_quartile3_bleach_peak_JulAug.csv"))

### Fisher's exact test: bleached vs unbleached across Q1-Q4 (July-August peak) ###
# 4 x 2 table per species (quartile x bleached / unbleached)
fisher_quartile_ja <- field_peak_ja %>%
  mutate(bleached = factor(bleach_max > 0, levels = c(TRUE, FALSE),
                           labels = c("Bleached", "Unbleached"))) %>%
  group_by(species) %>%
  group_modify(~ {
    tab <- table(.x$quartile, .x$bleached)
    tibble(n          = sum(tab),
           n_bleached = sum(tab[, "Bleached"]),
           bleached_Q1 = paste0(tab["Q1", "Bleached"], "/", sum(tab["Q1", ])),
           bleached_Q2 = paste0(tab["Q2", "Bleached"], "/", sum(tab["Q2", ])),
           bleached_Q3 = paste0(tab["Q3", "Bleached"], "/", sum(tab["Q3", ])),
           bleached_Q4 = paste0(tab["Q4", "Bleached"], "/", sum(tab["Q4", ])),
           p          = fisher.test(tab)$p.value)
  }) %>%
  ungroup()
fisher_quartile_ja
write_csv(fisher_quartile_ja, here("output", "tables", "fisher_quartile_bleach_peak_JulAug.csv"))

### Combined: Jul-Aug peak stacked bars (top) + ED50 group ranges (bottom) ###
fig_peak_ja_combined <-
  (fig_peak_ja_stack + labs(x = NULL)) /
  fig_group3_range +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(panel.spacing.x = unit(1.5, "lines"))
fig_peak_ja_combined
ggsave(here("output", "figures", "fig_quartile3_combined_peak_JulAug.png"), fig_peak_ja_combined,
       width = 9, height = 7, dpi = 300)
ggsave(here("output", "figures", "fig_quartile3_combined_peak_JulAug.pdf"), fig_peak_ja_combined,
       width = 9, height = 7)

# -----------------------------------------------------------------------------
# 7. Seasonal bleaching index (BI), May–November
# -----------------------------------------------------------------------------
# Source: field monitoring data (all tagged colonies, not only those with an
# ED50 value).
# BI = (0*C6 + 1*C5 + 2*C4 + 3*C3 + 4*C2 + 5*C1) / 5
# Ck = proportion of surveyed colonies with colour score k (1-6).
# BI ranges from 0 (all colour 6, healthy) to 1 (all colour 1, fully bleached).
fe <- field_raw %>%
  select(month = Month, location = Location, species = Species,
         colony = Colony, colour = Colour) %>%
  filter(!is.na(colour)) %>%
  mutate(month = factor(month, levels = month_order))

# BI from colour-score counts, within the grouping variables given in `...`
calc_bi <- function(data, ...) {
  bi <- data %>%
    count(..., colour) %>%
    group_by(...) %>%
    mutate(n_colony = sum(n), prop = n / n_colony) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = colour, values_from = prop,
                names_prefix = "C", names_sort = TRUE, values_fill = 0)
  # add any colour class absent from the data so every C1-C6 column exists
  for (k in paste0("C", 1:6)) if (!k %in% names(bi)) bi[[k]] <- 0
  bi %>%
    mutate(BI = (0 * C6 + 1 * C5 + 2 * C4 + 3 * C3 + 4 * C2 + 5 * C1) / 5) %>%
    select(..., n_colony, C1:C6, BI) %>%
    arrange(...)
}

### BI per species and month (all sites pooled) ###
bi_species <- calc_bi(fe, species, month)
print(bi_species, n = Inf)
write_csv(bi_species, here("output", "tables", "bleaching_index_species.csv"))
bi_species
### BI per species, site and month ###
# cells with only one colony surveyed are omitted (BI from a single colony)
bi_site <- calc_bi(fe, species, location, month)

bi_site %>% filter(n_colony == 1) %>% select(species, location, month, n_colony, BI)

bi_site <- bi_site %>% filter(n_colony > 1)
print(bi_site, n = Inf)
write_csv(bi_site, here("output", "tables", "bleaching_index_site.csv"))

### Plot: species BI per month, mean ± SE across sites ###
month_labels <- c("5" = "May", "6" = "Jun", "7" = "Jul", "8" = "Aug",
                  "9" = "Sep", "10" = "Oct", "11" = "Nov")
bi_species_cols <- c("G. pectinata" = "#0072B2", "P. speciosa" = "#F0AB00FF",
                     "P. acuta" = "#CC79A7", "P. lutea" = "purple")  # as in Part5

bi_species_site <- bi_site %>%
  group_by(species, month) %>%
  summarise(n_site = n(), n_colony = sum(n_colony),
            BI_mean = mean(BI), BI_se = sd(BI) / sqrt(n_site),
            .groups = "drop") %>%
  mutate(SpeciesLabel = factor(species_labels[species],
                               levels = names(bi_species_cols)),
         MonthNum = match(as.character(month), month.name))
print(bi_species_site, n = Inf)
write_csv(select(bi_species_site, -SpeciesLabel, -MonthNum),
          here("output", "tables", "bleaching_index_species_site_se.csv"))

dodge <- position_dodge(width = 0.3)
fig_bi_site <- ggplot(bi_species_site,
                      aes(MonthNum, BI_mean, colour = SpeciesLabel)) +
  geom_line(position = dodge, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(BI_mean - BI_se, 0), ymax = BI_mean + BI_se),
                width = 0.2, position = dodge) +
  geom_point(size = 2.5, position = dodge) +
  scale_colour_manual(values = bi_species_cols) +
  scale_x_continuous(breaks = 5:11, labels = month_labels) +
  labs(x = "Month", y = "Bleaching index (mean ± SE across sites)",
       colour = "Species") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))
fig_bi_site
ggsave(here("output", "figures", "fig_bleaching_index_site_se.png"), fig_bi_site,
       width = 7, height = 4.5, dpi = 300)
ggsave(here("output", "figures", "fig_bleaching_index_site_se.pdf"), fig_bi_site,
       width = 7, height = 4.5)
