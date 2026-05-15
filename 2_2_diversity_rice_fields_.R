# =============================================================================
# Script:  2_2_diversity_rice_fields_.R
# Purpose: Compare macroinvertebrate density, species richness, Shannon index,
#          and gamma diversity within rice paddy fields across three dataset
#          subsets: merged (center + ditch), center only, and ditch-only fields.
#          Includes abundance- and sample-based rarefaction curves (iNEXT),
#          alpha diversity plots, and LMMs/GLMMs for richness, Shannon index,
#          and density. Produces diversity figures for the manuscript.
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
library(ggh4x)
library(vegan)
library(iNEXT)
library(GGally)
library(wesanderson)
library(lme4)
library(lmerTest)
library(car)
library(sjPlot)
library(broom.mixed)
library(purrr)
library(DHARMa)

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

# --- 3.1 Add ditch presence column -------------------------------------------

ditch_yes <- c("brugg", "detligen", "jonen", "muehlau", "untersiggenthal_ost", "vionnaz")
ditch_no  <- c("la_sauge", "stetten", "untersiggenthal_west", "witzwil", "wuerenlingen")

comm_data1 <- comm_data %>%
  mutate(
    ditch_present = case_when(
      site_name %in% ditch_yes ~ "yes",
      site_name %in% ditch_no  ~ "no",
      TRUE                     ~ NA_character_
    ),
    .after = site_name
  )

# Subsampling area (quadrat: 25 × 25 cm)
subsample_area_m2 <- 0.25 * 0.25

# Filter to rice fields only and pivot to long format
comm_data_long <- comm_data1 %>%
  filter(env == "rice") %>%
  pivot_longer(
    cols      = Coenangrionidae:Sialis,
    names_to  = "species",
    values_to = "count"
  )


# --- 3.2 Full dataset: center and ditch kept separate (sp_data_f) ------------

sp_data_f <- comm_data_long %>%
  group_by(site_name, site_in, year, sampling_round, species) %>%
  summarise(
    count                   = mean(count / total_subamples) / subsample_area_m2,
    total_subsamples        = sum(total_subamples),
    ditch_present           = first(ditch_present),
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

meta_cols_f <- c(
  "site_name", "site_in", "year", "sampling_round",
  "total_subsamples", "ditch_present", "canton", "ch_NS",
  "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy",
  "ph", "temp", "NO", "PO", "N_fertiliser_name",
  "fertiliser_type", "total_N_per_ha", "fertiliser_latest_date",
  "soil_prep_type", "soil_prep_depth_cm", "soil_prep_date",
  "rice_transplant", "rice_harvest", "water_overwinter",
  "water_start", "water_stop"
)

# Remove species columns that are all zeros or NA
# Use dplyr::select explicitly to avoid conflicts with other packages
sp_data_f_cleaned <- sp_data_f %>%
  dplyr::select(all_of(meta_cols_f), where(~ any(. != 0, na.rm = TRUE)))

sp_data_fu <- sp_data_f_cleaned[, 30:ncol(sp_data_f_cleaned)]
env_f       <- sp_data_f_cleaned[, 1:29]

# Log(x + 1) transformation for area, conductivity, dissolved oxygen,
# water depth, temperature, and distance to wetland
env_f[, c(9:13, 15, 23)] <- log1p(env_f[, c(9:13, 15, 23)])

# Z-score scaling for ordination use only — do NOT use for LMMs or richness
env_f[, c(9:15, 23)] <- scale(env_f[, c(9:15, 23)])

env_f$N_fertiliser_name <- as.factor(env_f$N_fertiliser_name)
env_f$fertiliser_type   <- as.factor(env_f$fertiliser_type)
env_f$soil_prep_type    <- as.factor(env_f$soil_prep_type)
env_f$rice_transplant   <- as.factor(env_f$rice_transplant)
env_f$rice_harvest      <- as.factor(env_f$rice_harvest)

env_only_f <- env_f[, c(
  "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy",
  "ph", "temp", "N_fertiliser_name", "fertiliser_type", "total_N_per_ha",
  "soil_prep_type", "soil_prep_depth_cm", "rice_transplant", "rice_harvest"
)]


# --- 3.3 Merged dataset: center and ditch averaged per site (sp_data_m1) -----

sp_data_m1 <- comm_data_long %>%
  group_by(site_name, year, sampling_round, species) %>%
  summarise(
    count                   = mean(count / total_subamples) / subsample_area_m2,
    total_subsamples        = sum(total_subamples),
    ditch_present           = first(ditch_present),
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

meta_cols_m <- c(
  "site_name", "year", "sampling_round",
  "total_subsamples", "ditch_present", "canton", "ch_NS",
  "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy",
  "ph", "temp", "NO", "PO", "N_fertiliser_name",
  "fertiliser_type", "total_N_per_ha", "fertiliser_latest_date",
  "soil_prep_type", "soil_prep_depth_cm", "soil_prep_date",
  "rice_transplant", "rice_harvest", "water_overwinter",
  "water_start", "water_stop"
)

sp_data_m_cleaned <- sp_data_m1 %>%
  dplyr::select(all_of(meta_cols_m), where(~ any(. != 0, na.rm = TRUE)))

sp_data_m <- sp_data_m_cleaned[, 29:ncol(sp_data_m_cleaned)]
env_m      <- sp_data_m_cleaned[, 1:28]

# Log(x + 1) and z-score scaling for ordination use only
env_m[, c(8, 10:12, 14)] <- log1p(env_m[, c(8, 10:12, 14)])
env_m[, c(8, 10:14)]     <- scale(env_m[, c(8, 10:14)], center = TRUE)

env_m$year           <- as.factor(env_m$year)
env_m$sampling_round <- as.factor(env_m$sampling_round)
env_m$ch_NS          <- as.factor(env_m$ch_NS)

env_m_only <- env_m[, c(8, 10:14)]


# --- 3.4 Center-only dataset (sp_data_r) -------------------------------------

sp_data_r1 <- sp_data_f %>% filter(site_in == "center")

sp_data_r         <- sp_data_r1[, 30:ncol(sp_data_r1)]
env_r             <- sp_data_r1[, 1:29]
env_r[, c(9, 11:13, 15)] <- log1p(env_r[, c(9, 11:13, 15)])
env_r[, c(9, 11:15)]     <- scale(env_r[, c(9, 11:15)])

env_only_r <- env_r[, c(9, 11:15)]

# Cleaned center-only metadata + species matrix (site_in dropped via meta_cols_m)
sp_data_r_cleaned <- sp_data_r1 %>%
  dplyr::select(all_of(meta_cols_m), where(~ any(. != 0, na.rm = TRUE)))


# --- 3.5 Ditch-only dataset: fields with a ditch (sp_data_ditch) -------------

sp_data_ditch <- comm_data_long %>%
  filter(ditch_present == "yes") %>%
  group_by(site_name, site_in, year, sampling_round, species) %>%
  summarise(
    count                   = mean(count / total_subamples) / subsample_area_m2,
    total_subsamples        = sum(total_subamples),
    ditch_present           = first(ditch_present),
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

meta_cols <- c(
  "site_name", "site_in", "year", "sampling_round",
  "total_subsamples", "ditch_present", "canton", "ch_NS",
  "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy",
  "ph", "temp", "NO", "PO", "N_fertiliser_name",
  "fertiliser_type", "total_N_per_ha", "fertiliser_latest_date",
  "soil_prep_type", "soil_prep_depth_cm", "soil_prep_date",
  "rice_transplant", "rice_harvest", "water_overwinter",
  "water_start", "water_stop"
)

sp_data_ditch_cleaned <- sp_data_ditch %>%
  dplyr::select(all_of(meta_cols), where(~ any(. != 0, na.rm = TRUE)))

sp_data_d <- sp_data_ditch_cleaned[, 30:ncol(sp_data_ditch_cleaned)]
env_d      <- sp_data_ditch_cleaned[, 1:29]
env_d[, c(9, 11:13, 15)] <- log1p(env_d[, c(9, 11:13, 15)])
env_d[, c(9, 11:15)]     <- scale(env_d[, c(9, 11:15)])

env_only_d <- env_d[, c(9, 11:15)]


# -----------------------------------------------------------------------------
# 4. Macroinvertebrate density, richness, and Shannon index
# -----------------------------------------------------------------------------

# --- 4.1 Density per replicate (site × position × year × sampling round) -----

rep_density <- comm_data_long %>%
  group_by(site_name, site_in, year, sampling_round, species) %>%
  mutate(density = (count / total_subamples) / subsample_area_m2) %>%
  group_by(site_name, site_in, year, sampling_round) %>%
  summarise(
    total_ind_m2  = sum(density),
    ditch_present = first(ditch_present),
    .groups = "drop"
  ) %>%
  mutate(sampling_period = paste0(year, "_", sampling_round))

# Average density across sampling periods per site
site_density <- rep_density %>%
  group_by(site_name, site_in) %>%
  summarise(
    mean_ind_m2   = mean(total_ind_m2),
    ditch_present = first(ditch_present),
    .groups = "drop"
  )


# --- 4.2 Richness and Shannon index per replicate ----------------------------

alpha_data <- sp_data_f[, c(
  "site_name", "ditch_present", "site_in", "year", "sampling_round",
  names(sp_data_f)[30:ncol(sp_data_f)]
)]

alpha_pointwise <- alpha_data %>%
  group_by(site_name, site_in, year, sampling_round) %>%
  summarise(
    richness      = specnumber(across(where(is.numeric))),
    shannon       = diversity(across(where(is.numeric)), index = "shannon"),
    ditch_present = first(ditch_present),
    .groups = "drop"
  )

# Average richness and Shannon index across sampling periods per site
alpha_site_avg <- alpha_pointwise %>%
  group_by(site_name, site_in) %>%
  summarise(
    mean_alpha    = mean(richness),
    sd_alpha      = sd(richness),
    n             = n(),
    se_alpha      = sd_alpha / sqrt(n),
    mean_shannon  = mean(shannon),
    sd_shannon    = sd(shannon),
    se_shannon    = sd_shannon / sqrt(n),
    ditch_present = first(ditch_present),
    .groups = "drop"
  )


# --- 4.3 Combined alpha + density plot ---------------------------------------

combined_data <- left_join(
  alpha_site_avg, site_density,
  by = c("site_name", "site_in", "ditch_present")
)

custom_colors <- c(
  "Present - Side"   = "goldenrod",
  "Present - Center" = "gold2",
  "Absent - Center"  = "gray38"
)

plot_data <- combined_data %>%
  dplyr::mutate(
    ditch_present = dplyr::recode(as.character(ditch_present), "yes" = "Present", "no" = "Absent"),
    fill_group = case_when(
      ditch_present == "Present" & site_in == "side"   ~ "Present - Side",
      ditch_present == "Present" & site_in == "center" ~ "Present - Center",
      ditch_present == "Absent"  & site_in == "center" ~ "Absent - Center",
      TRUE ~ NA_character_
    )
  )

vars     <- c("mean_ind_m2", "mean_alpha", "mean_shannon")
y_labels <- c(
  "Macroinvertebrate \ndensity (ind/m²)",
  "Macroinvertebrate richness",
  "Shannon index"
)

alpha_plot_list <- list()

for (i in seq_along(vars)) {
  v    <- vars[i]
  ylab <- y_labels[i]

  alpha_plot_list[[v]] <- ggplot(
    plot_data, aes(x = site_in, y = .data[[v]], fill = fill_group, colour = fill_group)
  ) +
    geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
    stat_summary(
      aes(shape = fill_group),
      fun = mean, geom = "point", size = 4
    ) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.1, linewidth = 0.6) +
    scale_fill_manual(values = custom_colors, na.value = "white") +
    scale_colour_manual(values = custom_colors, na.value = "white") +
    scale_shape_manual(values = c(
      "Present - Side"   = 17,
      "Present - Center" = 19,
      "Absent - Center"  = 21
    )) +
    labs(
      x      = NULL,
      y      = ylab,
      fill   = "Ditch presence\n& sampling position",
      colour = "Ditch presence\n& sampling position",
      shape  = "Ditch presence\n& sampling position"
    ) +
    facet_grid2(
      . ~ ditch_present,
      scales = "free_x",
      space  = "free",
      strip  = strip_themed(
        background_x = elem_list_rect(fill = c("grey78", "palegoldenrod"))
      )
    ) +
    theme_pubr()
}

final_alpha_plot <- wrap_plots(alpha_plot_list, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top", legend.box = "vertical")

final_alpha_plot


# -----------------------------------------------------------------------------
# 5. Gamma diversity — rarefaction and extrapolation curves (iNEXT)
# -----------------------------------------------------------------------------

# Summed counts per site × year × sampling round (used for iNEXT)
sp_data_aggr <- comm_data_long %>%
  group_by(site_name, year, sampling_round, ditch_present, species) %>%
  summarise(
    count            = sum(count),
    total_subsamples = sum(total_subamples),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = "species", values_from = "count")

sp_data_aggr_center <- comm_data_long %>%
  filter(site_in == "center") %>%
  group_by(site_name, year, sampling_round, ditch_present, species) %>%
  summarise(
    count            = sum(count),
    total_subsamples = sum(total_subamples),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = "species", values_from = "count")

sp_data_aggr_ditch <- comm_data_long %>%
  filter(ditch_present == "yes") %>%
  group_by(site_name, site_in, year, sampling_round, species) %>%
  summarise(
    count            = sum(count),
    total_subsamples = sum(total_subamples),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = "species", values_from = "count")


# --- 5.1 Abundance-based rarefaction — merged dataset ------------------------

gamma_data <- sp_data_aggr[, c("site_name", "ditch_present", names(sp_data_aggr)[6:ncol(sp_data_aggr)])]

gamma_long <- gamma_data %>%
  group_by(ditch_present) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  column_to_rownames("ditch_present")

gamma_list <- apply(gamma_long, 1, function(x) x[x > 0])
gamma_list <- lapply(gamma_list, unlist)

out_gamma <- iNEXT(gamma_list, q = 0, datatype = "abundance")

p_gamma1 <- ggiNEXT(out_gamma, type = 1) +
  scale_shape_manual(
    name   = "Ditch presence in paddy",
    values = c(no = 19, yes = 15),
    labels = c(no = "Absent", yes = "Present")
  ) +
  scale_fill_manual(
    name   = "Ditch presence in paddy",
    values = c(no = "grey60", yes = "lightgoldenrod"),
    labels = c(no = "Absent", yes = "Present")
  ) +
  scale_colour_manual(
    name   = "Ditch presence in paddy",
    values = c(no = "grey38", yes = "lightgoldenrod3"),
    labels = c(no = "Absent", yes = "Present")
  ) +
  labs(y = "Macroinvertebrate\nrichness", tag = "A") +
  theme_pubr()


# --- 5.2 Abundance-based rarefaction — center only ---------------------------

gamma_data_C <- sp_data_aggr_center[, c("site_name", "ditch_present", names(sp_data_aggr_center)[6:ncol(sp_data_aggr_center)])]

gamma_long_C <- gamma_data_C %>%
  group_by(ditch_present) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  column_to_rownames("ditch_present")

gamma_list_C <- apply(gamma_long_C, 1, function(x) x[x > 0])
gamma_list_C <- lapply(gamma_list_C, unlist)

out_gamma_C <- iNEXT(gamma_list_C, q = 0, datatype = "abundance")

p_gamma1_C <- ggiNEXT(out_gamma_C, type = 1) +
  scale_shape_manual(
    name   = "Ditch presence in paddy",
    values = c(no = 15, yes = 16),
    labels = c(no = "Absent", yes = "Present")
  ) +
  scale_fill_manual(
    name   = "Ditch presence in paddy",
    values = c(no = "grey38", yes = "goldenrod"),
    labels = c(no = "Absent", yes = "Present")
  ) +
  scale_colour_manual(
    name   = "Ditch presence in paddy",
    values = c(no = "grey38", yes = "goldenrod"),
    labels = c(no = "Absent", yes = "Present")
  ) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  labs(y = "Macroinvertebrate\nrichness", tag = "A") +
  theme_pubr()


# --- 5.3 Abundance-based rarefaction — ditch-only fields ---------------------

gamma_data_D <- sp_data_aggr_ditch[, c("site_name", "site_in", names(sp_data_aggr_ditch)[6:ncol(sp_data_aggr_ditch)])]

gamma_long_D <- gamma_data_D %>%
  group_by(site_in) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  column_to_rownames("site_in")

gamma_list_D <- apply(gamma_long_D, 1, function(x) x[x > 0])
gamma_list_D <- lapply(gamma_list_D, unlist)

out_gamma_D <- iNEXT(gamma_list_D, q = 0, datatype = "abundance")

p_gamma1_D <- ggiNEXT(out_gamma_D, type = 1) +
  scale_shape_manual(
    name   = "Sampling position",
    values = c(center = 19, side = 17),
    labels = c(center = "Center", side = "Ditch")
  ) +
  scale_fill_manual(
    name   = "Sampling position",
    values = c(center = "lightgoldenrod", side = "burlywood3"),
    labels = c(center = "Center", side = "Ditch")
  ) +
  scale_colour_manual(
    name   = "Sampling position",
    values = c(center = "gold2", side = "goldenrod3"),
    labels = c(center = "Center", side = "Ditch")
  ) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  labs(y = "Macroinvertebrate\nrichness", tag = "B") +
  theme_pubr()

# Combined gamma panel (merged + ditch-only)
p_new_gamma <- p_gamma1 + p_gamma1_D +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

p_new_gamma / final_alpha_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none")


# --- 5.4 Sample-based rarefaction — merged dataset ---------------------------

comm_matrix <- sp_data_aggr[, -c(1:5)]

d_present        <- subset(comm_matrix, sp_data_aggr$ditch_present == "yes")
d_present2       <- d_present
d_present2[d_present2 > 0] <- 1
id_present       <- as.data.frame(as.incfreq(t(d_present2)))

d_absent         <- subset(comm_matrix, sp_data_aggr$ditch_present == "no")
d_absent2        <- d_absent
d_absent2[d_absent2 > 0] <- 1
id_absent        <- as.data.frame(as.incfreq(t(d_absent2)))

iditch           <- cbind(id_present, id_absent)
colnames(iditch) <- c("yes", "no")

ditch0 <- iNEXT(iditch, datatype = "incidence_freq", q = 0)

p_gamma2 <- ggiNEXT(ditch0, type = 1) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(no = "grey38", yes = "goldenrod")) +
  scale_colour_manual(values = c(no = "grey38", yes = "goldenrod")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  labs(y = "Macroinvertebrate diversity", x = "Number of samples", tag = "B") +
  theme_pubr()

p_gamma22 <- ggiNEXT(ditch0, type = 2) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(no = "grey38", yes = "goldenrod")) +
  scale_colour_manual(values = c(no = "grey38", yes = "goldenrod")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  theme_pubr()

p_gamma23 <- ggiNEXT(ditch0, type = 3) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(no = "grey38", yes = "goldenrod")) +
  scale_colour_manual(values = c(no = "grey38", yes = "goldenrod")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  theme_pubr()

p_gamma2 + p_gamma22 + p_gamma23


# --- 5.5 Sample-based rarefaction — center only ------------------------------

comm_matrix_C <- sp_data_aggr_center[, -c(1:5)]

d_present_C        <- subset(comm_matrix_C, sp_data_aggr_center$ditch_present == "yes")
d_present_C2       <- d_present_C
d_present_C2[d_present_C2 > 0] <- 1
id_present_C       <- as.data.frame(as.incfreq(t(d_present_C2)))

d_absent_C         <- subset(comm_matrix_C, sp_data_aggr_center$ditch_present == "no")
d_absent_C2        <- d_absent_C
d_absent_C2[d_absent_C2 > 0] <- 1
id_absent_C        <- as.data.frame(as.incfreq(t(d_absent_C2)))

iditch_C           <- cbind(id_present_C, id_absent_C)
colnames(iditch_C) <- c("yes", "no")

ditch0_C <- iNEXT(iditch_C, datatype = "incidence_freq", q = 0)

p_gamma2_C <- ggiNEXT(ditch0_C, type = 1) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(no = "grey38", yes = "goldenrod")) +
  scale_colour_manual(values = c(no = "grey38", yes = "goldenrod")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  labs(y = "Macroinvertebrate diversity", x = "Number of samples") +
  theme_pubr()

p_gamma2_C2 <- ggiNEXT(ditch0_C, type = 2) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(no = "grey38", yes = "goldenrod")) +
  scale_colour_manual(values = c(no = "grey38", yes = "goldenrod")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  theme_pubr()

p_gamma2_C3 <- ggiNEXT(ditch0_C, type = 3) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(no = "grey38", yes = "goldenrod")) +
  scale_colour_manual(values = c(no = "grey38", yes = "goldenrod")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  theme_pubr()

p_gamma2_C + p_gamma2_C2 + p_gamma2_C3 +
  plot_layout(guides = "collect") &
  ggtitle("Dataset: Center only") &
  theme(legend.position = "top")


# --- 5.6 Sample-based rarefaction — ditch-only fields ------------------------

comm_matrix_D <- sp_data_aggr_ditch[, -c(1:5)]

d_side        <- subset(comm_matrix_D, sp_data_aggr_ditch$site_in == "side")
d_side2       <- d_side
d_side2[d_side2 > 0] <- 1
id_side       <- as.data.frame(as.incfreq(t(d_side2)))

d_center      <- subset(comm_matrix_D, sp_data_aggr_ditch$site_in == "center")
d_center2     <- d_center
d_center2[d_center2 > 0] <- 1
id_center     <- as.data.frame(as.incfreq(t(d_center2)))

iditch_D           <- cbind(id_side, id_center)
colnames(iditch_D) <- c("side", "center")

ditch0_D <- iNEXT(iditch_D, datatype = "incidence_freq", q = 0)

p_gamma2_D <- ggiNEXT(ditch0_D, type = 1) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(center = "grey38", side = "green4")) +
  scale_colour_manual(values = c(center = "grey38", side = "green4")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  labs(x = "Number of samples", tag = "B") +
  theme_pubr()

p_gamma2_D2 <- ggiNEXT(ditch0_D, type = 2) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(center = "grey38", side = "green4")) +
  scale_colour_manual(values = c(center = "grey38", side = "green4")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  theme_pubr()

p_gamma2_D3 <- ggiNEXT(ditch0_D, type = 3) +
  scale_shape_manual(values = c(15, 16)) +
  scale_fill_manual(values   = c(center = "grey38", side = "green4")) +
  scale_colour_manual(values = c(center = "grey38", side = "green4")) +
  scale_y_continuous(breaks = seq(0, 60, 20)) +
  theme_pubr()

p_gamma2_D + p_gamma2_D2 + p_gamma2_D3 +
  plot_layout(guides = "collect") &
  ggtitle("Dataset: Ditch fields only") &
  theme(legend.position = "top")


# --- 5.7 Combined gamma panels -----------------------------------------------

gamma_panel   <- p_gamma1 + p_gamma2 &
  ggtitle("Dataset: Merged (center + ditch)") &
  theme(legend.position = "top")

gamma_panel_C <- p_gamma1_C + p_gamma2_C &
  ggtitle("Dataset: Center only") &
  theme(legend.position = "top")

gamma_panel_D <- p_gamma1_D + p_gamma2_D +
  plot_layout(guides = "collect") &
  ggtitle("Dataset: Ditch fields only") &
  theme(legend.position = "top")

gamma_m_c_panel <- gamma_panel / gamma_panel_C +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

gamma_m_c_panel / gamma_panel_D &
  plot_annotation(tag_levels = "A")

# Figure 5: combined gamma + alpha panel
gamma_panel / final_alpha_plot +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")


# -----------------------------------------------------------------------------
# 6. Descriptive statistics and pairwise scatterplots
# -----------------------------------------------------------------------------

env_data <- sp_data_f[, c(
  "distance_rice_wetland_m", "ch_NS", "area_ha", "oxy", "cond", "depth",
  "ph", "temp", "rice_transplant", "soil_prep_depth_cm",
  "N_fertiliser_name", "fertiliser_type",
  "site_name", "year", "sampling_round", "ditch_present", "site_in"
)]

env_data$transplant_doy <- as.numeric(format(env_data$rice_transplant, "%j"))
env_data <- env_data[, !(names(env_data) %in% "rice_transplant")]

species_data  <- sp_data_f_cleaned[, 30:ncol(sp_data_f_cleaned)]
combined_data <- cbind(env_data, species_data)
complete_rows <- complete.cases(combined_data)

env_clean     <- env_data[complete_rows, ]
env_clean_raw <- env_data[complete_rows, ]
species_clean <- species_data[complete_rows, ]

env_clean$sampling_period <- interaction(env_clean$year, env_clean$sampling_round)

# Per-group summary statistics for environmental predictors
summary_fields    <- list()
vars_to_summarize <- c(
  "distance_rice_wetland_m", "ch_NS", "area_ha", "oxy", "cond", "depth",
  "ph", "temp", "transplant_doy", "soil_prep_depth_cm",
  "N_fertiliser_name", "fertiliser_type"
)

for (v in vars_to_summarize) {
  var_value <- env_clean[[v]]

  if (is.numeric(var_value)) {
    summary_table <- env_clean %>%
      group_by(ditch_present, site_in, sampling_period) %>%
      summarise(
        variable = v,
        n        = sum(!is.na(.data[[v]])),
        mean     = mean(.data[[v]], na.rm = TRUE),
        sd       = sd(.data[[v]], na.rm = TRUE),
        .groups  = "drop"
      )

  } else if (is.character(var_value) || is.factor(var_value)) {
    summary_table <- env_clean %>%
      group_by(ditch_present, site_in, sampling_period, value = .data[[v]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(variable = v)

  } else if (inherits(var_value, "Date") || inherits(var_value, "POSIXct")) {
    summary_table <- env_clean %>%
      group_by(ditch_present, site_in, sampling_period) %>%
      summarise(
        variable = v,
        min_date = min(.data[[v]], na.rm = TRUE),
        max_date = max(.data[[v]], na.rm = TRUE),
        .groups  = "drop"
      )

  } else {
    next
  }

  summary_fields[[v]] <- summary_table
}

# Bray-Curtis dissimilarity (square-root transformed)
bray_dist <- sqrt(vegdist(species_clean, method = "bray"))

# Pairwise scatterplots for selected abiotic predictors, coloured by ditch presence
vars_fert   <- c("oxy", "cond", "ph", "temp", "N_fertiliser_name", "fertiliser_type")
vars_to_plot <- c(
  "distance_rice_wetland_m", "area_ha", "oxy", "cond", "depth", "ph", "temp",
  "transplant_doy", "soil_prep_depth_cm", "N_fertiliser_name", "fertiliser_type"
)

wes_col <- wes_palette("FantasticFox1", type = "discrete")

env_scatter <- ggpairs(
  env_clean_raw,
  columns = vars_fert,
  mapping = aes(colour = ditch_present, fill = ditch_present, alpha = 0.5),
  upper   = list(continuous = wrap(ggally_cor, size = 3)),
  lower   = list(continuous = wrap("smooth", method = "lm", se = TRUE)),
  diag    = list(continuous = "densityDiag"),
  legend  = 1
) +
  scale_fill_manual(values  = c(wes_col[1], wes_col[3])) +
  scale_color_manual(values = c(wes_col[1], wes_col[3])) +
  theme_pubr()

env_scatter + theme(axis.text.x = element_text(angle = 90, hjust = 1.0, vjust = 0.5))

ggpairs(
  env_clean_raw,
  columns = vars_to_plot,
  upper   = list(continuous = wrap(ggally_cor, size = 3)),
  lower   = list(continuous = wrap("smooth", method = "lm", se = TRUE)),
  diag    = list(continuous = "densityDiag")
) +
  theme_pubr()


# =============================================================================
# 7. Linear Mixed Models (LMMs / GLMMs)
# =============================================================================

# Three dataset subsets are modelled:
#   r: center only,         fixed effect = ditch_present
#   m: merged (center+ditch) averaged, fixed effect = ditch_present
#   d: ditch-only fields,   fixed effect = site_in (center vs. ditch)
#
# For r and m datasets, ditch_present is a between-site factor, so sites are
# nested within ditch categories: (1 | ditch_present:site_name).
# For d, site_in is a within-site factor: (1 | site_name).

# Named list of datasets with their cleaned metadata and species start column
datasets <- list(
  r = list(data = sp_data_r,     meta = sp_data_r_cleaned,     start_col = 30),
  m = list(data = sp_data_m,     meta = sp_data_m_cleaned,     start_col = 29),
  d = list(data = sp_data_ditch, meta = sp_data_ditch_cleaned, start_col = 30)
)

plot_titles <- list(
  r = "Center only — all fields",
  m = "Merged (center + ditch) — all fields",
  d = "Ditch-only fields"
)

# Type III ANOVA requires sum-to-zero contrasts for correct hypothesis tests
# with unbalanced designs
options(contrasts = c("contr.treatment", "contr.poly"))


combined_results <- list()
model_results <- list()
anova_results <- list()
model_plots <- list()


for (name in names(datasets)) {
  
  sp_data_cleaned <- datasets[[name]]$meta
  start_col <- datasets[[name]]$start_col
  species_cols <- colnames(sp_data_cleaned)[start_col:ncol(sp_data_cleaned)]
  
  grouping_vars <- c("site_name", "year", "sampling_round")
  if ("site_in" %in% colnames(sp_data_cleaned)) {
    grouping_vars <- c(grouping_vars, "site_in")
  }
  
  alpha_data <- sp_data_cleaned %>%
    dplyr::select(all_of(grouping_vars), ditch_present, all_of(species_cols)) %>%
    mutate(sampling_period = interaction(year, sampling_round, drop = TRUE))
  
  result <- alpha_data %>%
    group_by(across(all_of(grouping_vars)), sampling_period) %>%
    reframe(
      richness = specnumber(as.matrix(across(all_of(species_cols)))),
      shannon = diversity(as.matrix(across(all_of(species_cols))), index = "shannon"),
      ditch_present = first(ditch_present),
      density = sum(across(all_of(species_cols)), na.rm = TRUE),
      dataset = name
    )
  
  combined_results[[name]] <- result
  
  # --- Fixed effect variable ---
  if (name %in% c("m", "r")) {
    fixed_effect <- "ditch_present"
  } else if (name == "d") {
    fixed_effect <- "site_in"
  } else {
    fixed_effect <- NULL
  }
  
  if (!is.null(fixed_effect)) {
    
    # --- Random effect structure ---
    if (name %in% c("m", "r")) {
      # ditch_present is between-site → site is nested inside ditch category
      rand_effects <- "(1|ditch_present:site_name) + (1|sampling_period)"
      # use the explicit one if you want site random intercepts preserved
    } else if (name == "d") {
      # site_in is within-site → keep site random intercept
      rand_effects <- "(1|site_name) + (1|sampling_period)"
    }
    
    # --- Models ---
    mod_richness <- tryCatch({
      glmer(
        as.formula(paste0("richness ~ factor(", fixed_effect, ") + ", rand_effects)),
        data = result,
        family = poisson(link = "log")
      )
    }, error = function(e) e)
    
    mod_shannon <- tryCatch({
      lmer(
        as.formula(paste0("shannon ~ factor(", fixed_effect, ") + ", rand_effects)),
        data = result
      )
    }, error = function(e) e)
    
    mod_density <- tryCatch({
      lmer(
        as.formula(paste0("log1p(density) ~ factor(", fixed_effect, ") + ", rand_effects)),
        data = result
      )
    }, error = function(e) e)
    
    model_results[[name]] <- list(
      richness_model = mod_richness,
      shannon_model = mod_shannon,
      density_model = mod_density
    )
    
    # --- ANOVA ---
    anova_results[[name]] <- list(
      richness_anova = if (name == "d") {
        Anova(mod_richness, type = "III") # GLMM
      } else {
        Anova(mod_richness, type = "III")
      },
      shannon_anova  = anova(mod_shannon, type = "III"),
      density_anova  = anova(mod_density, type = "III")
    )
    
    # --- Plots ---
    model_plots[[name]] <- plot_models(
      model_results[[name]]$richness_model,
      model_results[[name]]$shannon_model,
      model_results[[name]]$density_model,
      show.values = TRUE,
      show.p = TRUE,
      title = plot_titles[[name]],
      legend.title = "Response",
      axis.labels = c(paste(fixed_effect, sep = "_"))
    )
  }
}




# -----------------------------------------------------------------------------
# 8. Model outputs
# -----------------------------------------------------------------------------

# --- 8.1 Model summary tables ------------------------------------------------

tab_model(model_results$m)
tab_model(model_results$r)
tab_model(model_results$d)


# --- 8.2 ANOVA tables --------------------------------------------------------

anova(model_results$d$density_model)
anova(model_results$d$shannon_model)
Anova(model_results$d$richness_model, type = 3)

anova(model_results$m$density_model)
anova(model_results$m$shannon_model)
Anova(model_results$m$richness_model, type = 3)


# Combined ANOVA table with BH-adjusted p-values (via broom.mixed)

tidy_anova_normalized <- function(anova_obj) {
  df <- tidy(anova_obj)
  
  if ("NumDF" %in% names(df)) {
    # LMM anova(): F-test with num/den df
    df %>%
      rename(num.df = NumDF, den.df = DenDF, f.value = statistic) %>%
      mutate(stat_type = "F", chisq = NA_real_, df_chisq = NA_integer_) %>%
      dplyr::select(term, stat_type, f.value, chisq, num.df, den.df, df_chisq, sumsq, p.value)
  } else {
    # GLMM Anova(): chi-square test with single df
    df %>%
      rename(chisq = statistic, df_chisq = df) %>%
      mutate(stat_type = "Chisq", f.value = NA_real_,
             num.df = NA_integer_, den.df = NA_real_, sumsq = NA_real_) %>%
      dplyr::select(term, stat_type, f.value, chisq, num.df, den.df, df_chisq, sumsq, p.value)
  }
}

anova_table <- purrr::imap_dfr(anova_results, function(dataset_models, dataset_name) {
  purrr::imap_dfr(dataset_models, function(anova_obj, response_name) {
    tidy_anova_normalized(anova_obj) %>%
      mutate(dataset = dataset_name, response = response_name)
  })
}) %>%
  group_by(dataset) %>%
  mutate(
    p.adj   = p.adjust(p.value, method = "BH"),
    p.value = format.pval(p.value, digits = 3, eps = 0.001),
    p.adj   = format.pval(p.adj,   digits = 3, eps = 0.001)
  ) %>%
  ungroup()

anova_table


# --- 8.3 DHARMa residual diagnostics -----------------------------------------

plot_dharma_dataset <- function(dataset_name, n_sims = 250) {
  mods <- model_results[[dataset_name]]
  if (is.null(mods)) stop("No models for dataset: ", dataset_name)

  graphics.off()
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

  make_panel <- function(mod, title_txt) {
    if (inherits(mod, "error")) {
      plot.new()
      title(main = paste(title_txt, "\n(model failed)"))
      return(invisible(NULL))
    }
    sim <- DHARMa::simulateResiduals(mod, n = n_sims)
    plot(sim, main = title_txt)
  }

  make_panel(mods$richness_model, paste(dataset_name, "- Richness"))
  make_panel(mods$shannon_model,  paste(dataset_name, "- Shannon"))
  make_panel(mods$density_model,  paste(dataset_name, "- Density"))
}

plot_dharma_dataset("r")
plot_dharma_dataset("m")
plot_dharma_dataset("d")


# --- 8.4 Coefficient plots ---------------------------------------------------

# Per-response plots across all three datasets
model_plots_response <- list(richness = list(), shannon = list(), density = list())

for (name in names(model_results)) {
  fixed_effect <- ifelse(name %in% c("m", "r"), "ditch_present", "site_in")

  model_plots_response$richness[[name]] <- plot_models(
    model_results[[name]]$richness_model,
    show.values  = TRUE, show.p = TRUE,
    title        = paste0(plot_titles[[name]], " — Richness (IRR)"),
    transform    = "exp",
    legend.title = "Response",
    axis.labels  = fixed_effect
  )

  model_plots_response$shannon[[name]] <- plot_models(
    model_results[[name]]$shannon_model,
    show.values  = TRUE, show.p = TRUE,
    title        = paste0(plot_titles[[name]], " — Shannon"),
    legend.title = "Response",
    axis.labels  = fixed_effect
  )

  model_plots_response$density[[name]] <- plot_models(
    model_results[[name]]$density_model,
    show.values  = TRUE, show.p = TRUE,
    title        = paste0(plot_titles[[name]], " — Density"),
    legend.title = "Response",
    axis.labels  = fixed_effect
  )
}

# Combined coefficient plots across datasets per response variable
richness_modplot <- plot_models(
  model_results$d$richness_model,
  model_results$r$richness_model,
  model_results$m$richness_model,
  title      = "Richness",
  m.labels   = c(
    "Position in field \n(Ditch only fields)",
    "Field center only",
    "Merged \n(Ditch + center)"
  ),
  show.values = TRUE,
  show.p      = TRUE,
  wrap.labels = 40
)

shannon_modplot <- plot_models(
  model_results$d$shannon_model,
  model_results$r$shannon_model,
  model_results$m$shannon_model,
  title      = "Shannon",
  m.labels   = c(
    "Position in field \n(Ditch only fields)",
    "Field center only",
    "Merged \n(Ditch + center)"
  ),
  show.values = TRUE,
  show.p      = TRUE,
  wrap.labels = 40
)

density_modplot <- plot_models(
  model_results$d$density_model,
  model_results$r$density_model,
  model_results$m$density_model,
  title      = "Log+1(Density)",
  m.labels   = c(
    "Position in field \n(Ditch only fields)",
    "Field center only",
    "Merged \n(Ditch + center)"
  ),
  show.values = TRUE,
  show.p      = TRUE,
  wrap.labels = 40
)

richness_modplot / shannon_modplot / density_modplot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top") &
  theme_pubr()

# Combined per-dataset coefficient plots
model_plots$r / model_plots$m / model_plots$d +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") &
  theme_pubr()
