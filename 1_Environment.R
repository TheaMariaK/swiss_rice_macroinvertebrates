# =============================================================================
# Script:  1_Environment.R
# Purpose: Compare abiotic environmental conditions between rice paddy fields
#          and reference wetlands. Produces descriptive statistics, assumption
#          tests, MANOVA, LDA, and LMMs for individual environmental variables.
#          Produces Figure 1 of the manuscript.
# Author:  Thea Bulas
# Manuscript title: Rice paddies promote diverse and distinct aquatic
#                   invertebrate communities in agroecosystems.
# Date:    15.05.2026
# Data:    2022_2023_macroinvRADICAL_abio_env_20250407.xlsx
# =============================================================================   


# -----------------------------------------------------------------------------
# 1. Libraries
# ----------------------------------------------------------------------------- 

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggdist)
library(patchwork)
library(readxl)
library(corrplot)
library(rstatix)
library(mvnormalTest)
library(heplots)
library(MASS)
library(effectsize)
library(lme4)
library(lmerTest)
library(gridExtra)
library(car)
library(purrr)
library(GGally)
library(lmerTest)
library(gridExtra)
library(car)
library(purrr)
library(GGally)
library(sjPlot)

setwd("~/RStudio/agroscope/chapter_2/Manuscript_scripts")


# -----------------------------------------------------------------------------
# 2. Data import
# -----------------------------------------------------------------------------

comm_data <- read_excel(
  "~/RStudio/agroscope/chapter_2/Manuscript_scripts/data_and_Rscripts_CLEAN/data/2022_2023_macroinvRADICAL_abio_env_20250407.xlsx",
  col_types = c("text", "text", "numeric", 
                "numeric", "numeric", "date", "text", 
                "numeric", "text", "numeric", "numeric", 
                "text", "text", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "text", "text", "numeric", 
                "date", "text", "date", "numeric", 
                "date", "date", "text", "date", "date", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric"))



# -----------------------------------------------------------------------------
# 3. Data wrangling
# -----------------------------------------------------------------------------

# Pivot species counts to long format
comm_data_long <- comm_data |>
  pivot_longer(
    cols      = Coenangrionidae:Sialis,
    names_to  = "species",
    values_to = "count"
  )

# Subsampling area (quadrat: 25 × 25 cm)
subsample_area_m2 <- 0.25 * 0.25

# Aggregate to site × year × sampling round × species level.
# Invertebrate density expressed as individuals per m².
# Abiotic variables averaged across within-site replicates.
sp_data_m1 <- comm_data_long |>
  group_by(site_name, env, year, sampling_round, species) |>
  summarise(
    count                   = mean(count / total_subamples) / subsample_area_m2,
    total_subsamples        = sum(total_subamples),
    canton                  = first(canton),
    ch_NS                   = first(ch_NS),
    area_ha                 = sum(area_ha),
    distance_rice_wetland_m = mean(distance_rice_wetland_m),
    # Abiotic measurements
    cond                    = mean(cond),
    depth                   = mean(depth),
    oxy                     = mean(oxy),
    ph                      = mean(ph),
    temp                    = mean(temp),
    NO                      = mean(NO),
    PO                      = mean(PO),
    # Fertiliser management
    N_fertiliser_name       = first(N_fertiliser_name),
    fertiliser_type         = first(fertiliser_type),
    total_N_per_ha          = mean(total_N_per_ha),
    fertiliser_latest_date  = min(fertiliser_latest_date),
    # Soil preparation
    soil_prep_type          = first(soil_prep_type),
    soil_prep_depth_cm      = mean(soil_prep_depth_cm),
    soil_prep_date          = min(soil_prep_date),
    # Rice phenology
    rice_transplant         = min(rice_transplant),
    rice_harvest            = min(rice_harvest),
    # Water management
    water_overwinter        = min(water_overwinter),
    water_start             = min(water_start),
    water_stop              = min(water_stop),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = "species", values_from = "count")


# -----------------------------------------------------------------------------
# 4. Environmental data preparation
# -----------------------------------------------------------------------------

# Subset to metadata and abiotic variables; convert grouping columns to factors
env_m <- sp_data_m1[, 1:16]

# Log(x + 1) transformation for conductivity, dissolved oxygen, water depth,
# temperature, and water body area
env_m[, c(8, 10:12, 14)] <- log1p(env_m[, c(8, 10:12, 14)])

env_m$env            <- as.factor(env_m$env)
env_m$year           <- as.factor(env_m$year)
env_m$sampling_round <- as.factor(env_m$sampling_round)
env_m$ch_NS          <- as.factor(env_m$ch_NS)

# Abiotic variables only (no metadata), used for multivariate analyses
env_only <- env_m[, c(8, 10:14)]


# -----------------------------------------------------------------------------
# 5. Descriptive statistics
# -----------------------------------------------------------------------------

# Per-habitat summary statistics (n, mean, SD) for each abiotic variable,
# computed on untransformed values
summary_list <- list()

for (v in colnames(sp_data_m1[, c(8, 10:14)])) {
  summary_table <- sp_data_m1 |>
    group_by(env = sp_data_m1$env) |>
    summarise(
      variable = v,
      n        = n(),
      mean     = mean(.data[[v]], na.rm = TRUE),
      sd       = sd(.data[[v]], na.rm = TRUE),
      .groups  = "drop"
    )
  summary_list[[v]] <- summary_table
}

summary_df <- do.call(rbind, summary_list)


# -----------------------------------------------------------------------------
# 6. Assumption testing
# -----------------------------------------------------------------------------

# --- 6.1 Univariate normality (Shapiro–Wilk test, per habitat) ---------------

result <- env_m |>
  group_by(env) |>
  shapiro_test(area_ha, cond, depth, oxy, ph, temp)

result |>
  mutate(p = format(p, scientific = FALSE, digits = 3)) |>
  print(n = Inf)


# --- 6.2 Multivariate normality (Mardia's test) -------------------------------

mardia(env_only[, c("area_ha", "cond", "depth", "oxy", "ph", "temp")])$mv.test


# --- 6.3 Homogeneity of covariance matrices (Box's M-test) -------------------

boxM(Y = env_only, group = env_m$env)


# --- 6.4 Levene's test for homogeneity of variance (per variable) ------------

vars <- names(env_only)

results_levene <- lapply(vars, function(v) {
  test <- car::leveneTest(as.formula(paste(v, "~ env")), data = env_m)
  data.frame(variable = v, statistic = test$`F value`[1], p = test$`Pr(>F)`[1])
})

levene_df <- do.call(rbind, results_levene)
levene_df


# --- 6.5 Multicollinearity (Pearson correlation matrix) ----------------------

cor(env_only, method = "pearson")


# --- 6.6 Linearity (pairwise scatterplots per habitat group) -----------------

env_plot_data <- sp_data_m1[, c(2, 8, 10:14)]  # includes env grouping column

ggpairs(
  env_plot_data,
  columns = 2:7,
  mapping = aes(colour = env, fill = env, alpha = 0.5),
  upper   = list(continuous = wrap(ggally_cor, size = 3)),
  lower   = list(continuous = wrap("smooth", method = "lm", se = TRUE)),
  diag    = list(continuous = "densityDiag")
) +
  theme_minimal()


# -----------------------------------------------------------------------------
# 7. Multivariate outlier detection and removal
# -----------------------------------------------------------------------------

# Flag outliers using Mahalanobis distance; remove before MANOVA and LDA
# One outlier identified: Jonen 2023, sampling round 2
outliers <- mahalanobis_distance(env_only)$is.outlier

env_m_no_outlier    <- env_m[!outliers, ]
env_only_no_outlier <- env_only[!outliers, ]


# -----------------------------------------------------------------------------
# 8. MANOVA
# -----------------------------------------------------------------------------

# Test the effect of habitat type on the multivariate set of abiotic variables
manova_fit <- manova(
  as.matrix(env_only_no_outlier) ~ env,
  data = env_m_no_outlier
)

summary(manova_fit, test = "Pillai")
summary(manova_fit, test = "Wilks")

# Univariate follow-up ANOVAs per response variable
manova.terms.aov <- summary.aov(manova_fit)

# Partial eta squared (effect size)
eta_squared(manova_fit)


# -----------------------------------------------------------------------------
# 9. Linear Discriminant Analysis (LDA)
# -----------------------------------------------------------------------------

# Identify the linear combination of abiotic variables that best discriminates
# between habitat types; extract and plot LD1 scores
lda_fit    <- lda(env ~ as.matrix(env_only_no_outlier), data = env_m_no_outlier)
lda_scores <- predict(lda_fit)$x
lda_plot   <- data.frame(env = env_m_no_outlier$env, lda_scores)


# --- 9.1 LDA plots -----------------------------------------------------------

# Boxplot of LD1 scores by habitat type
p_lda_overall <- ggplot(lda_plot, aes(x = env, y = LD1, fill = env, colour = env)) +
  geom_boxplot(width = 0.3, alpha = 0.2) +
  geom_jitter(width = 0.1) +
  scale_colour_manual(name = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_fill_manual(name   = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  labs(title = "A") +
  theme_pubr() +
  theme(axis.title.x = element_blank())

# Density plot of LD1 scores (overall)
lda_dens_overall <- ggplot(lda_plot, aes(x = LD1, fill = env, colour = env)) +
  geom_density(alpha = 0.2, show.legend = FALSE) +
  scale_colour_manual(name = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_fill_manual(name   = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  labs(title = "B") +
  theme_pubr()

# Density plot of LD1 scores faceted by year × sampling round
# Note: row order of lda_plot matches env_m_no_outlier
lda_dens_time <- ggplot(lda_plot, aes(x = LD1, fill = env, colour = env)) +
  geom_density(alpha = 0.2) +
  scale_colour_manual(values = c("#d55e00", "#0072b2")) +
  scale_fill_manual(values   = c("#d55e00", "#0072b2")) +
  labs(title = "LD1 by sampling period") +
  facet_wrap(~ env_m_no_outlier$year * env_m_no_outlier$sampling_round) +
  theme_pubr()

# Combined LDA panel (boxplot + overall density)
p_lda_overall + lda_dens_overall +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")


# -----------------------------------------------------------------------------
# 10. Linear Mixed Models (LMMs) for individual abiotic variables
# -----------------------------------------------------------------------------

# Habitat type (env) as fixed effect; site identity and sampling period
# (year × sampling round) as crossed random intercepts
env_m_no_outlier$sampling_period <- interaction(
  env_m_no_outlier$year,
  env_m_no_outlier$sampling_round
)

# LMMs for water-column variables (measured at each sampling round)
predictors <- c("cond", "depth", "oxy", "ph", "temp")
model_list <- list()

for (pred in predictors) {
  form               <- as.formula(paste(pred, "~ env + (1 | site_name) + (1 | sampling_period)"))
  model_list[[pred]] <- lmer(form, data = env_m_no_outlier, REML = TRUE)
}

# Separate LMM for water body area (time-invariant; site random effect only)
mod_area <- lmer(area_ha ~ env + (1 | site_name), data = env_m_no_outlier, REML = TRUE)

# Type III ANOVA tables (Satterthwaite denominator df)
anova(mod_area)
anova(model_list[["cond"]])
anova(model_list[["depth"]])
anova(model_list[["oxy"]])
anova(model_list[["ph"]])
anova(model_list[["temp"]])

# Model summary table (fixed effects, SE, confidence intervals)
tab_model(
  mod_area,
  model_list[["cond"]],
  model_list[["depth"]],
  model_list[["oxy"]],
  model_list[["ph"]],
  model_list[["temp"]],
  show.ngroups = FALSE,
  show.icc     = FALSE,
  show.zeroinf = FALSE,
  show.obs     = FALSE,
  show.se      = TRUE,
  collapse.ci  = TRUE
)

# Coefficient plot across all LMMs
plot_models(
  mod_area,
  model_list[["cond"]],
  model_list[["depth"]],
  model_list[["oxy"]],
  model_list[["ph"]],
  model_list[["temp"]],
  show.intercept = TRUE
) +
  theme_pubr()


# -----------------------------------------------------------------------------
# 11. Figure 1 — Abiotic conditions by habitat type
# -----------------------------------------------------------------------------

# Significance labels from LMM ANOVA results; y.position at 110% of observed max
stat_tbl <- tibble(
  variable   = c("area_ha", "cond", "depth", "oxy", "ph", "temp"),
  label      = c("p < 0.001", "p < 0.001", "n.s.", "n.s.", "n.s.", "p = 0.037"),
  group1     = "rice",
  group2     = "wetland",
  y.position = map_dbl(variable, ~ max(sp_data_m1[[.x]], na.rm = TRUE) * 1.1)
)

vars   <- c("area_ha", "cond", "depth", "oxy", "ph", "temp")
titles <- rep("", length(vars))
ytag   <- c(
  "Area (ha)", "Conductivity (µS/cm)", "Depth (cm)",
  "Dissolved oxygen (mg/L)", "pH", "Temperature (°C)"
)

# Build one panel per variable (raincloud plot: half-violin + boxplot)
plots <- pmap(
  list(vars, titles, ytag),
  function(v, title, ylabel) {
    ggplot(sp_data_m1, aes(env, .data[[v]], fill = env, colour = env)) +
      stat_halfeye(
        adjust        = 0.6,
        width         = 0.6,
        .width        = 0,
        justification = 0,
        point_colour  = NA,
        alpha         = 0.25,
        show.legend   = FALSE
      ) +
      geom_boxplot(
        width         = 0.20,
        outlier.shape = NA,
        alpha         = 0.25,
        show.legend   = FALSE
      ) +
      scale_colour_manual(name = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
      scale_fill_manual(name   = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
      labs(title = title, x = NULL, y = ylabel) +
      stat_pvalue_manual(
        data          = filter(stat_tbl, variable == v),
        mapping       = aes(x = group1, xend = group2, y = y.position),
        label         = "label",
        tip.length    = 0,
        step.increase = 0,
        inherit.aes   = FALSE,
        size          = 3.5
      ) +
      coord_cartesian(clip = "off") +
      theme_pubr() +
      theme(
        axis.title.y = element_text(size = 11),
        plot.margin  = margin(t = 15, r = 5, b = 5, l = 5)
      )
  }
)

# Figure 1: six-panel raincloud plot (panels A–F), shared legend at top
combined_plot <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")

combined_plot