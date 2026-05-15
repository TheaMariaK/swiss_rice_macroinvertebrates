# =============================================================================
# Script:  2_diversity.R
# Purpose: Estimate and compare gamma and alpha diversity (species richness,
#          Shannon index) and macroinvertebrate density between rice paddy
#          fields and reference wetlands. Includes abundance- and
#          sample-based rarefaction curves (iNEXT), LMMs for richness and
#          Shannon index, and a GLMM for species richness. Produces the
#          diversity figure for the manuscript.
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
library(patchwork)
library(readxl)
library(vegan)
library(iNEXT)
library(lme4)
library(lmerTest)
library(DHARMa)
library(sjPlot)
library(car)
library(multilevelTools)
library(JWileymisc)

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
comm_data_long <- comm_data %>%
  pivot_longer(
    cols      = Coenangrionidae:Sialis,
    names_to  = "species",
    values_to = "count"
  )

# Subsampling area (quadrat: 25 × 25 cm)
subsample_area_m2 <- 0.25 * 0.25

# Summed counts per site × habitat × year × sampling round (used for iNEXT)
sp_data_aggr <- comm_data_long %>%
  group_by(site_name, env, year, sampling_round, species) %>%
  summarise(
    count                   = sum(count),
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
  ) %>%
  pivot_wider(names_from = "species", values_from = "count")

# Density-based aggregation (ind/m²) averaged across subsamples
sp_data_m1 <- comm_data_long %>%
  group_by(site_name, env, year, sampling_round, species) %>%
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
  ) %>%
  pivot_wider(names_from = "species", values_from = "count")


# -----------------------------------------------------------------------------
# 4. Gamma diversity — rarefaction and extrapolation curves (iNEXT)
# -----------------------------------------------------------------------------

# --- 4.1 Abundance-based rarefaction (Hill number q = 0) ---------------------

# Aggregate total counts per habitat across all sites and sampling periods
gamma_data <- sp_data_aggr[, c("site_name", "env", names(sp_data_aggr)[29:ncol(sp_data_aggr)])]

gamma_long <- gamma_data %>%
  group_by(env) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  mutate(across(where(is.numeric), ceiling)) %>%
  column_to_rownames("env")

# Convert to named list for iNEXT input
gamma_list <- apply(gamma_long, 1, function(x) x[x > 0])
gamma_list <- lapply(gamma_list, unlist)

# Run abundance-based rarefaction and extrapolation
out_gamma <- iNEXT(gamma_list, q = 0, datatype = "abundance")

p_gamma1 <- ggiNEXT(out_gamma, type = 1) +
  scale_shape_manual(
    name   = "Habitat",
    values = c(rice = 16, wetland = 15),
    labels = c(rice = "Rice", wetland = "Wetland")
  ) +
  scale_fill_manual(
    name   = "Habitat",
    values = c(rice = "#d55e00", wetland = "#0072b2"),
    labels = c(rice = "Rice", wetland = "Wetland")
  ) +
  scale_colour_manual(
    name   = "Habitat",
    values = c(rice = "#d55e00", wetland = "#0072b2"),
    labels = c(rice = "Rice", wetland = "Wetland")
  ) +
  scale_y_continuous(breaks = seq(0, 60, 20), limits = c(0, 60)) +
  labs(y = "Macroinvertebrate diversity", tag = "A") +
  guides(
    shape  = guide_legend(title = "Habitat"),
    fill   = guide_legend(title = "Habitat"),
    colour = guide_legend(title = "Habitat")
  ) +
  theme_pubr()


# --- 4.2 Sample-based rarefaction (incidence frequencies, q = 0) -------------

# Prepare incidence frequency matrices per habitat
comm_matrix <- sp_data_m1[, -c(1:28)]

rice        <- subset(comm_matrix, sp_data_m1$env == "rice")
rice2       <- rice
rice2[rice2 > 0] <- 1
rice2       <- t(rice2)
irice       <- as.data.frame(as.incfreq(rice2))

wetland     <- subset(comm_matrix, sp_data_m1$env == "wetland")
wetland2    <- wetland
wetland2[wetland2 > 0] <- 1
wetland2    <- t(wetland2)
iwetland    <- as.data.frame(as.incfreq(wetland2))

ihab              <- cbind(irice, iwetland)
colnames(ihab)    <- c("rice", "wetland")

# Run sample-based rarefaction and extrapolation
hab0 <- iNEXT(ihab, datatype = "incidence_freq", q = 0)

p_gamma2 <- ggiNEXT(hab0, type = 1) +
  scale_shape_manual(
    name   = "Habitat",
    values = c(rice = 16, wetland = 15),
    labels = c(rice = "Rice", wetland = "Wetland")
  ) +
  scale_fill_manual(
    name   = "Habitat",
    values = c(rice = "#d55e00", wetland = "#0072b2"),
    labels = c(rice = "Rice", wetland = "Wetland")
  ) +
  scale_colour_manual(
    name   = "Habitat",
    values = c(rice = "#d55e00", wetland = "#0072b2"),
    labels = c(rice = "Rice", wetland = "Wetland")
  ) +
  scale_y_continuous(breaks = seq(0, 60, 20), limits = c(0, 60)) +
  labs(y = "Macroinvertebrate diversity", x = "Number of samples") +
  guides(
    shape  = guide_legend(title = "Habitat"),
    fill   = guide_legend(title = "Habitat"),
    colour = guide_legend(title = "Habitat")
  ) +
  theme_pubr()

gamma_panel <- p_gamma1 + p_gamma2


# -----------------------------------------------------------------------------
# 5. Alpha diversity — species richness and Shannon index
# -----------------------------------------------------------------------------

# Select species columns and add sampling period label
alpha_data <- sp_data_m1[, c("site_name", "env", "year", "sampling_round",
                              names(sp_data_m1)[29:ncol(sp_data_m1)])] %>%
  mutate(sampling_period = paste0(year, "_", sampling_round))

# Calculate richness and Shannon index per site × sampling period
alpha_pointwise <- alpha_data %>%
  group_by(site_name, env, sampling_period, year, sampling_round) %>%
  summarise(
    richness = specnumber(across(where(is.numeric))),
    shannon  = diversity(across(where(is.numeric)), index = "shannon"),
    .groups  = "drop"
  )

# Average richness and Shannon index per site across all sampling periods
alpha_site_avg <- alpha_pointwise %>%
  group_by(site_name, env) %>%
  summarise(
    mean_alpha   = mean(richness),
    sd_alpha     = sd(richness),
    n            = n(),
    se_alpha     = sd_alpha / sqrt(n),
    mean_shannon = mean(shannon),
    sd_shannon   = sd(shannon),
    se_shannon   = sd_shannon / sqrt(n),
    .groups = "drop"
  )


# --- 5.1 Richness plot (per site, averaged across sampling periods) ----------

alpha_sites <- ggplot(alpha_site_avg, aes(x = env, y = mean_alpha, fill = env, colour = env)) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  stat_summary(
    aes(colour = env, shape = env, fill = env),
    fun = mean, geom = "point", size = 5, show.legend = FALSE
  ) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.1, linewidth = 1, show.legend = FALSE) +
  labs(x = "", y = "Macroinvertebrate \nrichness", tag = "D") +
  scale_fill_manual(name   = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_colour_manual(name = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_shape_manual(values = c(rice = 21, wetland = 22)) +
  theme_pubr()


# --- 5.2 Shannon index plot --------------------------------------------------

shannon_sites <- ggplot(alpha_site_avg, aes(x = env, y = mean_shannon, fill = env, colour = env)) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  stat_summary(
    aes(colour = env, shape = env, fill = env),
    fun = mean, geom = "point", size = 5, show.legend = FALSE
  ) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.1, linewidth = 1, show.legend = FALSE) +
  labs(x = "", y = "Shannon index", tag = "E") +
  scale_fill_manual(name   = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_colour_manual(name = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_shape_manual(values = c(rice = 21, wetland = 22)) +
  theme_pubr()


# -----------------------------------------------------------------------------
# 6. Macroinvertebrate density
# -----------------------------------------------------------------------------

# Total density per replicate (site × year × sampling round)
rep_density <- comm_data_long %>%
  group_by(site_name, env, year, sampling_round) %>%
  mutate(density = sum(count / total_subamples / subsample_area_m2)) %>%
  summarise(mean_ind_m2 = mean(density), .groups = "drop") %>%
  mutate(sampling_period = paste0(year, "_", sampling_round))

# Average density across sampling periods per site
site_density <- rep_density %>%
  group_by(site_name, env) %>%
  summarise(mean_ind_m2 = mean(mean_ind_m2), .groups = "drop")

p_density <- ggplot(site_density, aes(x = env, y = mean_ind_m2, fill = env, colour = env)) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  stat_summary(
    aes(colour = env, shape = env, fill = env),
    fun = mean, geom = "point", size = 5, show.legend = FALSE
  ) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.1, linewidth = 1, show.legend = FALSE) +
  labs(x = NULL, y = "Macroinvertebrate \ndensity (ind/m\u00B2)", tag = "C") +
  scale_fill_manual(name   = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_colour_manual(name = "Habitat", values = c(rice = "#d55e00", wetland = "#0072b2")) +
  scale_shape_manual(values = c(rice = 21, wetland = 22)) +
  theme_pubr()


# -----------------------------------------------------------------------------
# 7. Combined figure
# -----------------------------------------------------------------------------

alpha_panel <- p_density + alpha_sites + shannon_sites

gamma_panel / alpha_panel +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")


# -----------------------------------------------------------------------------
# 8. Linear Mixed Models (LMMs / GLMMs)
# -----------------------------------------------------------------------------

# --- 8.1 Species richness (GLMM, Poisson) ------------------------------------

mod_richness <- glmer(
  richness ~ factor(env) + (1 | site_name) + (1 | sampling_period),
  data   = alpha_pointwise,
  family = poisson(link = "log")
)

Anova(mod_richness, type = "III")
tab_model(mod_richness)

# Residual diagnostics
simulationOutput_rich <- simulateResiduals(mod_richness)
plot(simulationOutput_rich)


# --- 8.2 Shannon index (LMM, Gaussian) ---------------------------------------

mod_shannon <- lmer(
  shannon ~ factor(env) + (1 | site_name) + (1 | sampling_period),
  data = alpha_pointwise
)

anova(mod_shannon)
tab_model(mod_shannon)

# Residual diagnostics
md <- modelDiagnostics(mod_shannon, ev.perc = .1)
plot(md, ncol = 2, nrow = 2)

simulationOutput_shannon <- simulateResiduals(mod_shannon)
plot(simulationOutput_shannon)


# --- 8.3 Macroinvertebrate density (LMM, log-transformed) -------------------

mod_dens <- lmer(
  log1p(mean_ind_m2) ~ factor(env) + (1 | site_name) + (1 | sampling_period),
  data = rep_density
)

anova(mod_dens)

# Residual diagnostics
simulationOutput_dens <- simulateResiduals(mod_dens)
plot(simulationOutput_dens)


# --- 8.4 Combined model outputs ----------------------------------------------

plot_models(mod_dens, mod_richness, mod_shannon) + theme_pubr()

tab_model(mod_dens, mod_richness, mod_shannon)
