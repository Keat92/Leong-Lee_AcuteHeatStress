# =============================================================================
# Part 2 — Acute heat-stress assay: dose–response curves and ED50
# =============================================================================
# Fits 4-parameter log-logistic dose–response curves (drc::LL.4) of
# photosynthetic efficiency (Fv/Fm, `y2`) against assay temperature for four
# coral species, then derives ED50 (temperature at 50 % loss of Fv/Fm):
#   * per species   — curve + confidence band for the figure,
#   * per fragment  — one curve per fragment (Tank x Colony),
#   * per colony    — mean of its fragment ED50s, ranked within species.
#
# Then tests whether fragment ED50 differs among species and sites with a
# linear mixed-effects model (colony random intercept, species-specific
# residual variances).
#
# Extracted from Part5_drc.Rmd (sections "Setup" to "Rank colonies within
# species", and "Statistical test" to the emmeans contrasts). The quartile
# analysis, variance partitioning and field comparison follow in later scripts.
#
# Inputs  (data/raw/):
#   consolidated.csv            PAM readings: Site, Time, Temp, Treatment,
#                               Tank, Colony, Species, F, Fm, Y(II)
# Outputs:
#   data/processed/fragment_ed50.csv       ED50 (+ SE) per fragment
#   data/processed/colony_ed50.csv         colony mean / SD / SE of ED50
#   data/processed/colony_ed50_ranked.csv  colony ED50 + within-species rank
#   output/tables/colony_ed50_range.csv    within-species colony ED50 spread
#   output/figures/fig2.tiff   species dose–response curves
#   output/figures/fig4.tiff   colony ED50 forest plot
#   output/figures/fig5.tiff   species ED50 box plot
#   output/figures/fig3.pdf    combined fig5 (A) + fig4 (B)
#   output/figures/ed50_model_diagnostics.pdf  assumption checks + model checks
#   output/tables/ed50_model_aic.csv           candidate variance structures
#   output/tables/ed50_anova_typeIII.csv       Type III Wald chi-square tests
#   output/tables/ed50_emmeans_site_by_species.csv  Tukey site contrasts
# =============================================================================

# NOTE: load `drc` BEFORE `dplyr`. `drc` attaches `MASS`, whose `select()`
# would otherwise mask `dplyr::select()` ("unused arguments" errors).
# Load order alone is not enough in a session where MASS was already attached
# (e.g. after running other scripts), so `select()` is also written as
# `dplyr::select()` below.
library(drc)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggsci)      # scale_colour_bmj(), scale_colour_aaas()
library(ggpubr)     # theme_pubr()
library(patchwork)  # combine panels
library(dlookr)     # describe()
library(nlme)       # lme() with varIdent() variance structures
library(car)        # leveneTest(), Anova()
library(performance)# check_model(), check_collinearity(), check_singularity()
library(see)        # plotting backend for check_model()
library(emmeans)    # post-hoc contrasts
library(here)

# -----------------------------------------------------------------------------
# 1. Load and format data
# -----------------------------------------------------------------------------
data_raw <- read.csv(here("data", "raw", "consolidated.csv"))
str(data_raw)

# Colonies re-measured in a second run carry a ".2" suffix; only the first
# measurement of each colony is analysed.
remeasured <- c("GOS09.2", "GOS10.2", "GOS11.2", "GOS12.2", "GOS13.2",
                "PAC11.2", "PAC12.2", "PAC14.2",
                "POC11.2", "POC12.2", "POC15.2",
                "POR11.2", "POR12.2", "POR13.2")

data <- data_raw %>%
  rename(y2 = X1.Y..II.,
         f  = X1.F,
         fm = X1.Fm) %>%
  mutate(
    Species   = as.factor(Species),
    Colony    = as.factor(Colony),
    Site      = as.factor(Site),
    Treatment = as.factor(Treatment),
    Fragment  = as.factor(paste(Tank, Colony, sep = "_"))
  ) %>%
  filter(Treatment == "Heated",
         !Colony %in% remeasured) %>%
  droplevels()

# Guard: no other suffixed colony IDs should survive the filter.
stopifnot(!any(grepl("\\.[0-9]+$", levels(data$Colony))))

# Sample sizes
data %>% summarise(n_fragments = n_distinct(Fragment),
                   n_colonies  = n_distinct(Colony))
data %>% group_by(Species) %>% summarise(n_colonies = n_distinct(Colony))
data %>% group_by(Site)    %>% summarise(n_colonies = n_distinct(Colony))

data %>% group_by(Species, Site, Colony) %>% describe()

# -----------------------------------------------------------------------------
# 2. Plot aesthetics
# -----------------------------------------------------------------------------
# Full binomial names for the species codes (facet strips and legends).
species_labels <- c(
  GOS = "G. pectinata",
  POC = "P. acuta",
  POR = "P. lutea",
  PAC = "P. speciosa"
)

species_cols <- c(
  "G. pectinata" = "#0072B2",
  "P. speciosa"  = "#F0AB00FF",
  "P. acuta"     = "#CC79A7",
  "P. lutea"     = "purple"
)

plot_family <- "sans"
plot_size   <- 12

# -----------------------------------------------------------------------------
# 3. Species-level dose–response curves
# -----------------------------------------------------------------------------
# One LL.4 model per species, used for the fitted curves in Fig. 2.
species_models <- lapply(split(data, data$Species), function(dat_sub) {
  drm(y2 ~ Temp, data = dat_sub,
      fct = LL.4(names = c("slope", "lower", "upper", "ED50")))
})

# Species ED50 from the pooled species curve (delta-method 95 % CI), for
# comparison with the fragment-based estimates below.
species_ed50_curve <- bind_rows(lapply(names(species_models), function(sp) {
  ed <- ED(species_models[[sp]], 50, interval = "delta", display = FALSE)
  tibble(Species = sp, ED50 = ed[1, "Estimate"],
         Lower = ed[1, "Lower"], Upper = ed[1, "Upper"])
}))
species_ed50_curve

temp_grid <- data.frame(
  Temp = seq(min(data$Temp), max(data$Temp), length.out = 200)
)

pred_df <- bind_rows(lapply(names(species_models), function(sp) {
  pr <- as.data.frame(
    predict(species_models[[sp]], newdata = temp_grid, interval = "confidence")
  )
  data.frame(Temp       = temp_grid$Temp,
             Prediction = pr[, 1],
             Lower      = pr[, 2],
             Upper      = pr[, 3],
             Species    = sp)
}))

# -----------------------------------------------------------------------------
# 4. Fragment-level ED50
# -----------------------------------------------------------------------------
# Each colony has >= 6 fragments (one per tank), so ED50 is estimated per
# fragment and then averaged to the colony.
fragment_model_list <- lapply(split(data, data$Fragment), function(dat_sub) {
  tryCatch(
    drm(y2 ~ Temp, data = dat_sub,
        fct = LL.4(names = c("slope", "lower", "upper", "ED50"))),
    error = function(e) NULL
  )
})

failed <- names(Filter(is.null, fragment_model_list))
message(length(failed), " of ", length(fragment_model_list),
        " fragment fits failed",
        if (length(failed)) paste0(": ", paste(failed, collapse = ", ")))

fragment_meta <- data %>%
  distinct(Fragment, Species, Site, Colony, Tank)

fragment_ed50 <- bind_rows(lapply(names(fragment_model_list), function(frag) {
  mod <- fragment_model_list[[frag]]
  if (is.null(mod)) return(NULL)
  cf <- summary(mod)$coefficients
  tibble(Fragment = frag,
         ED50     = cf["ED50:(Intercept)", "Estimate"],
         SE       = cf["ED50:(Intercept)", "Std. Error"])
})) %>%
  left_join(fragment_meta %>% mutate(Fragment = as.character(Fragment)),
            by = "Fragment") %>%
  dplyr::select(Fragment, Species, Site, Colony, Tank, ED50, SE)

fragment_ed50

write_csv(fragment_ed50, here("data", "processed", "fragment_ed50.csv"))

# -----------------------------------------------------------------------------
# 5. Colony-level ED50
# -----------------------------------------------------------------------------
colony_ed50 <- fragment_ed50 %>%
  group_by(Site, Species, Colony) %>%
  summarise(
    n_frag    = n(),
    ED50_mean = mean(ED50, na.rm = TRUE),
    ED50_sd   = sd(ED50, na.rm = TRUE),
    ED50_se   = ED50_sd / sqrt(n_frag),
    .groups   = "drop"
  )

write_csv(colony_ed50, here("data", "processed", "colony_ed50.csv"))

# Species mean ED50 (dashed lines in Fig. 2). Averaged over colony means so
# each colony counts once, regardless of how many fragments fitted.
ed50_sp_mean <- colony_ed50 %>%
  group_by(Species) %>%
  summarise(ED50_mean = mean(ED50_mean),
            ED50_sd   = sd(ED50_mean),
            n_colony  = n(),
            .groups   = "drop")
ed50_sp_mean

# Spread of colony-mean ED50 within each species
colony_ed50_range <- colony_ed50 %>%
  group_by(Species) %>%
  summarise(
    n_colony   = n(),
    ED50_min   = min(ED50_mean, na.rm = TRUE),
    ED50_max   = max(ED50_mean, na.rm = TRUE),
    ED50_range = ED50_max - ED50_min,
    ED50_IQR   = IQR(ED50_mean, na.rm = TRUE),
    Colony_min = Colony[which.min(ED50_mean)],
    Colony_max = Colony[which.max(ED50_mean)],
    .groups    = "drop"
  )
colony_ed50_range

write_csv(colony_ed50_range, here("output", "tables", "colony_ed50_range.csv"))

# Five-number summary per species (matches the box plot below)
species_ed50_summary <- colony_ed50 %>%
  group_by(Species) %>%
  summarise(
    ED50_min    = min(ED50_mean, na.rm = TRUE),
    ED50_q1     = quantile(ED50_mean, 0.25, na.rm = TRUE),
    ED50_median = median(ED50_mean, na.rm = TRUE),
    ED50_q3     = quantile(ED50_mean, 0.75, na.rm = TRUE),
    ED50_max    = max(ED50_mean, na.rm = TRUE),
    .groups     = "drop"
  )
species_ed50_summary

# Rank colonies within species: 1 = highest ED50 (most heat tolerant);
# ties share the lowest rank.
colony_ed50_ranked <- colony_ed50 %>%
  group_by(Species) %>%
  mutate(Rank = min_rank(desc(ED50_mean))) %>%
  arrange(Species, Rank) %>%
  ungroup()
colony_ed50_ranked

write_csv(colony_ed50_ranked, here("data", "processed", "colony_ed50_ranked.csv"))

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------
# Fig. 2 — species dose–response curves with mean colony ED50
p1 <- ggplot(data, aes(x = Temp, y = y2, color = Species)) +
  geom_point(alpha = 0.2) +
  geom_ribbon(
    data = pred_df,
    aes(x = Temp, ymin = Lower, ymax = Upper, fill = Species, group = Species),
    alpha = 0.2, inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_df,
    aes(x = Temp, y = Prediction, color = Species, group = Species),
    inherit.aes = FALSE, linewidth = 0.5
  ) +
  geom_vline(
    data = ed50_sp_mean,
    aes(xintercept = ED50_mean, color = Species),
    linetype = "dashed", show.legend = FALSE
  ) +
  geom_text(
    data = ed50_sp_mean,
    aes(x = ED50_mean, y = Inf, color = Species,
        label = paste0("ED50 = ", round(ED50_mean, 2))),
    angle = 90, vjust = -0.5, hjust = 1.1, inherit.aes = FALSE,
    show.legend = FALSE,
    family = plot_family, size = plot_size / .pt
  ) +
  scale_color_bmj(labels = species_labels) +
  scale_fill_bmj(labels = species_labels) +
  labs(x = "Temperature (°C)", y = "Fv/Fm",
       color = "Species", fill = "Species") +
  guides(fill = "none") +
  theme_pubr(base_size = plot_size, base_family = plot_family) +
  theme(legend.text = element_text(face = "italic")) +
  coord_cartesian(clip = "off")
p1

ggsave(here("output", "figures", "fig2.tiff"), plot = p1,
       width = 11.7, height = 8.3, dpi = 300, compression = "lzw")

# Fig. 4 — colony mean ED50 ± SD, ordered within species
p2 <- colony_ed50 %>%
  mutate(Colony = reorder(Colony, ED50_mean)) %>%
  ggplot(aes(x = ED50_mean, y = Colony, colour = Site)) +
  geom_errorbar(aes(xmin = ED50_mean - ED50_sd,
                    xmax = ED50_mean + ED50_sd),
                width = 0.2, alpha = 0.7, orientation = "y") +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ Species, scales = "free_y",
             labeller = as_labeller(species_labels)) +
  scale_colour_aaas() +
  labs(x = "Mean ED50 (°C) ± SD", y = "Colony") +
  theme_pubr(base_size = plot_size, base_family = plot_family) +
  theme(axis.text.y = element_text(size = 9),
        strip.text  = element_text(face = "italic"))
p2

ggsave(here("output", "figures", "fig4.tiff"), plot = p2,
       width = 11.7, height = 8.3, dpi = 300, compression = "lzw")

# Fig. 5 — species distribution of colony ED50. Box = IQR, line = median,
# whiskers (coef = Inf) = full min–max; species ordered by median.
p_range <- colony_ed50 %>%
  mutate(SpeciesLabel = reorder(species_labels[as.character(Species)],
                                ED50_mean, FUN = median)) %>%
  ggplot(aes(x = ED50_mean, y = SpeciesLabel, colour = SpeciesLabel)) +
  geom_boxplot(outlier.shape = NA, coef = Inf, width = 0.6, alpha = 0.5) +
  geom_jitter(height = 0.15, width = 0, size = 1.5, alpha = 0.5) +
  scale_colour_manual(values = species_cols) +
  guides(colour = "none") +
  labs(x = "ED50 (°C)", y = "Species") +
  theme_pubr(base_size = plot_size, base_family = plot_family) +
  theme(axis.text.y = element_text(face = "italic"))
p_range

ggsave(here("output", "figures", "fig5.tiff"), plot = p_range,
       width = 8, height = 5, dpi = 300, compression = "lzw")

# Fig. 3 — species range (A) above colony forest plot (B)
fig3 <- p_range / p2 +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "A")
fig3

ggsave(here("output", "figures", "fig3.pdf"), plot = fig3,
       width = 8.3, height = 11)

# -----------------------------------------------------------------------------
# 7. Statistical analysis of fragment ED50
# -----------------------------------------------------------------------------
# Response: fragment ED50 (6 fragments per colony). Fixed: Species * Site.
# Random: colony intercept (fragments of one colony are not independent).
stats_dat <- fragment_ed50 %>%
  mutate(across(c(Species, Site, Colony), as.factor))

# Every Species x Site cell must be present for the interaction to be estimable.
with(distinct(stats_dat, Species, Site, Colony), table(Species, Site))

## 7a. Assumption checks on raw ED50 ------------------------------------------
# Tests on the pooled raw response are a first look only; the residual checks
# on the fitted model (7c) are the ones that matter.
shapiro.test(stats_dat$ED50)

leveneTest(ED50 ~ Species,        data = stats_dat)
leveneTest(ED50 ~ Site,           data = stats_dat)
leveneTest(ED50 ~ Species * Site, data = stats_dat)

# Transformations: over a 33–38 °C range log() and sqrt() are almost linear,
# so they cannot be expected to fix skew or unequal variance.
stats_dat <- stats_dat %>%
  mutate(log_ED50  = log(ED50),
         sqrt_ED50 = sqrt(ED50))

shapiro.test(stats_dat$log_ED50)
shapiro.test(stats_dat$sqrt_ED50)
leveneTest(log_ED50  ~ Species * Site, data = stats_dat)
leveneTest(sqrt_ED50 ~ Species * Site, data = stats_dat)
# No improvement -> model raw ED50 and allow unequal variances instead.

## 7b. Linear mixed-effects models with unequal variances --------------------
# Same fixed and random parts; only the residual-variance structure differs,
# so the REML fits can be compared directly by AIC.
mod_none <- lme(ED50 ~ Species * Site, random = ~ 1 | Colony,
                data = stats_dat, method = "REML")
mod_sp   <- update(mod_none, weights = varIdent(form = ~ 1 | Species))
mod_si   <- update(mod_none, weights = varIdent(form = ~ 1 | Site))

model_aic <- anova(mod_none, mod_sp, mod_si)
model_aic
write_csv(as_tibble(model_aic, rownames = "model") %>%
            dplyr::select(-call),
          here("output", "tables", "ed50_model_aic.csv"))

# Selected: lowest AIC. Refit with ML for the fixed-effect tests.
# Sum-to-zero contrasts are required for meaningful Type III tests; they are
# set only now (REML logLik/AIC above depend on the contrast coding, so the
# candidate comparison stays on the default treatment contrasts) and restored
# at the end of this section. update() refits, so the new coding is used.
old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))
best_ML <- update(mod_sp, method = "ML")
summary(best_ML)

## 7c. Model diagnostics ------------------------------------------------------
check_collinearity(best_ML)
check_singularity(best_ML)

pdf(here("output", "figures", "ed50_model_diagnostics.pdf"),
    width = 10, height = 10)
hist(stats_dat$ED50, breaks = 30, main = "Fragment ED50", xlab = "ED50 (°C)")
qqnorm(stats_dat$ED50, main = "Normal Q-Q: raw ED50"); qqline(stats_dat$ED50)
print(check_model(best_ML))
dev.off()

## 7d. Fixed-effect tests and post-hoc contrasts ------------------------------
ed50_anova <- car::Anova(best_ML, type = "III")
ed50_anova
write_csv(as_tibble(as.data.frame(ed50_anova), rownames = "term"),
          here("output", "tables", "ed50_anova_typeIII.csv"))

# Site differences within each species (Tukey-adjusted)
ed50_emm <- emmeans(best_ML, pairwise ~ Site | Species, adjust = "tukey")
ed50_emm
write_csv(as_tibble(summary(ed50_emm$contrasts)),
          here("output", "tables", "ed50_emmeans_site_by_species.csv"))

options(old_contrasts)

sessionInfo()
