# =============================================================================
# Script:  1_2_Environment_ricefields.R
# Purpose: Compare abiotic environmental conditions within rice paddy fields
#          across three subsets: merged (center + ditch), center only, and
#          fields with a ditch only. Includes descriptive statistics, LMMs
#          for individual environmental variables. Produces
#          figures for the manuscript.
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
library(ggdist)
library(ggh4x)
library(lme4)
library(lmerTest)
library(sjPlot)
library(purrr)
library(MVN)
library(car)
library(rstatix)
library(effectsize)
library(MASS)
library(heplots)

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

# Aggregate to site × position × year × sampling round × species level.
# Density expressed as ind/m²; abiotic variables averaged across replicates.
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

# Separate species matrix and environmental metadata
sp_data_fu <- sp_data_f_cleaned[, 30:ncol(sp_data_f_cleaned)]
env_f       <- sp_data_f_cleaned[, 1:29]

# Log(x + 1) transformation for area, conductivity, dissolved oxygen,
# water depth, temperature, and distance to wetland (for ordinations only)
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

# Log(x + 1) transformation for conductivity, dissolved oxygen, water depth,
# temperature, and water body area
env_m[, c(8, 10:12, 14)] <- log1p(env_m[, c(8, 10:12, 14)])

env_m$site_name      <- as.factor(env_m$site_name)
env_m$ditch_present  <- as.factor(env_m$ditch_present)
env_m$year           <- as.factor(env_m$year)
env_m$sampling_round <- as.factor(env_m$sampling_round)
env_m$ch_NS          <- as.factor(env_m$ch_NS)

env_m_only <- env_m[, c(8, 10:14)]


# --- 3.4 Center-only dataset (sp_data_r) -------------------------------------

sp_data_r1 <- sp_data_f %>% filter(site_in == "center")

sp_data_r <- sp_data_r1[, 30:ncol(sp_data_r1)]
env_r      <- sp_data_r1[, 1:29]

env_r[, c(9, 11:13, 15)] <- log1p(env_r[, c(9, 11:13, 15)])

env_r$site_name      <- as.factor(env_r$site_name)
env_r$ditch_present  <- as.factor(env_r$ditch_present)
env_r$year           <- as.factor(env_r$year)
env_r$sampling_round <- as.factor(env_r$sampling_round)
env_r$ch_NS          <- as.factor(env_r$ch_NS)

env_only_r <- env_r[, c(9, 11:15)]


# --- 3.5 Ditch-only dataset: fields with a ditch, center and ditch (sp_data_ditch) ---

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

env_d$site_name      <- as.factor(env_d$site_name)
env_d$ditch_present  <- as.factor(env_d$ditch_present)
env_d$year           <- as.factor(env_d$year)
env_d$sampling_round <- as.factor(env_d$sampling_round)
env_d$ch_NS          <- as.factor(env_d$ch_NS)
env_d$site_in        <- as.factor(env_d$site_in)

env_only_d <- env_d[, c(9, 11:15)]


# -----------------------------------------------------------------------------
# 4. Descriptive statistics
# -----------------------------------------------------------------------------

# Per-habitat summary statistics (n, mean, SD) grouped by ditch presence and
# sampling position, for numeric variables; frequency counts for categorical;
# date range for date variables
summary_fields    <- list()
vars_to_summarize <- c("area_ha", "cond", "depth", "oxy", "ph", "temp")

for (v in vars_to_summarize) {
  var_value <- sp_data_f_cleaned[[v]]

  if (is.numeric(var_value)) {
    summary_table <- sp_data_f_cleaned %>%
      group_by(ditch_present, site_in) %>%
      summarise(
        variable = v,
        n        = sum(!is.na(.data[[v]])),
        mean     = mean(.data[[v]], na.rm = TRUE),
        sd       = sd(.data[[v]], na.rm = TRUE),
        .groups  = "drop"
      )

  } else if (is.character(var_value) || is.factor(var_value)) {
    summary_table <- sp_data_f_cleaned %>%
      group_by(ditch_present, site_in, value = .data[[v]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(variable = v)

  } else if (inherits(var_value, "Date") || inherits(var_value, "POSIXct")) {
    summary_table <- sp_data_f_cleaned %>%
      group_by(ditch_present, site_in) %>%
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

# Combined summary table for numeric variables
summary_df <- do.call(
  rbind,
  summary_fields[vars_to_summarize[sapply(sp_data_f[vars_to_summarize], is.numeric)]]
)
summary_df


# -----------------------------------------------------------------------------
# 5. Figures — environmental variables by ditch presence
# -----------------------------------------------------------------------------

# --- 5.1 Overall distribution (all sampling rounds combined) -----------------

# Recode ditch presence labels and assign fill groups for plotting
sp_data_f_cleaned <- sp_data_f_cleaned %>%
  mutate(
    ditch_present = dplyr::recode(ditch_present, "yes" = "Present", "no" = "Absent"),
    fill_group = case_when(
      ditch_present == "Present" & site_in == "side"   ~ "Present - Side",
      ditch_present == "Present" & site_in == "center" ~ "Present - Center",
      ditch_present == "Absent"  & site_in == "center" ~ "Absent - Center",
      TRUE ~ NA_character_
    )
  )

custom_colors <- c(
  "Present - Side"   = "#E69F00",
  "Present - Center" = "#F0E442",
  "Absent - Center"  = "gray40"
)

vars   <- c("area_ha", "cond", "depth", "oxy", "ph", "temp")
titles <- rep("", length(vars))
ytag   <- c(
  "Area (ha)", "Conductivity (µS/cm)", "Depth (cm)",
  "Dissolved oxygen (mg/L)", "pH", "Temperature (°C)"
)

plot_list <- list()

for (i in seq_along(vars)) {
  v <- vars[i]

  plot_list[[v]] <- ggplot(sp_data_f_cleaned, aes(x = site_in, y = .data[[v]], fill = fill_group)) +
    geom_boxplot(width = 0.20, outlier.shape = NA, alpha = 0.25) +
    stat_halfeye(adjust = 0.6, width = 0.6, .width = 0,
                 justification = 0, point_colour = NA, alpha = 0.25) +
    labs(
      title = titles[i],
      x     = "",
      y     = ytag[i],
      fill  = "Ditch presence - Measurement position:"
    ) +
    scale_fill_manual(values = custom_colors, na.value = "white") +
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

# Six-panel plot (panels A–F), shared legend at top
wrap_plots(plot_list, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")


# --- 5.2 Distribution over time (faceted by sampling period) -----------------

# Recode sampling labels and create ordered sampling period factor
sp_data_f_cleaned <- sp_data_f_cleaned %>%
  dplyr::mutate(
    ditch_present  = factor(
      dplyr::recode(as.character(ditch_present), "yes" = "Present", "no" = "Absent"),
      levels = c("Absent", "Present")
    ),
    sampling_round = dplyr::recode(as.character(sampling_round), "1" = "July",   "2" = "August"),
    year           = dplyr::recode(as.character(year),           "2022" = "'22", "2023" = "'23"),
    sampling_period = factor(
      paste(sampling_round, year),
      levels = c("July '22", "August '22", "July '23", "August '23")
    ),
    fill_group = case_when(
      ditch_present == "Present" & site_in == "side"   ~ "Present - Side",
      ditch_present == "Present" & site_in == "center" ~ "Present - Center",
      ditch_present == "Absent"  & site_in == "center" ~ "Absent - Center",
      TRUE ~ NA_character_
    )
  )

vars_t <- c("cond", "depth", "oxy", "ph", "temp")
ytag_t <- c(
  "Conductivity (µS/cm)", "Depth (cm)",
  "Dissolved oxygen (mg/L)", "pH", "Temperature (°C)"
)

plot_list_t <- list()

for (i in seq_along(vars_t)) {
  v <- vars_t[i]

  plot_list_t[[v]] <- ggplot(sp_data_f_cleaned, aes(x = site_in, y = .data[[v]], fill = fill_group)) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25) +
    stat_halfeye(adjust = 0.6, width = 0.6, .width = 0,
                 justification = 0, point_colour = NA, alpha = 0.25) +
    labs(
      x    = "",
      y    = ytag_t[i],
      fill = "Ditch presence - Measurement position:"
    ) +
    scale_fill_manual(values = custom_colors, na.value = "white") +
    facet_nested(
      cols       = vars(sampling_period, ditch_present),
      rows       = NULL,
      scales     = "free_x",
      space      = "free_x",
      nest_line  = TRUE,
      strip      = strip_nested(bleed = FALSE)
    ) +
    theme_pubr()
}

# Five-panel plot (panels A–E), shared legend at top
wrap_plots(plot_list_t, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position  = "top",
    strip.text.x     = element_text(margin = margin(0.1, 0, 0.1, 0, "cm"))
  )


# -----------------------------------------------------------------------------
# 6. Environmental summaries
# -----------------------------------------------------------------------------

summary(env_m)
summary(env_r)
summary(env_d)


# -----------------------------------------------------------------------------
# 7. Linear Mixed Models (LMMs) for individual abiotic variables
# -----------------------------------------------------------------------------

# Three dataset subsets are modelled:
#   env_m: merged (center + ditch averaged), fixed effect = ditch_present
#   env_r: center only,                      fixed effect = ditch_present
#   env_d: ditch-only fields,                fixed effect = site_in (center vs. ditch)
# Site identity and sampling round are included as crossed random intercepts.

datasets <- list(
  env_m = list(data = env_m, fixed = "ditch_present"),
  env_r = list(data = env_r, fixed = "ditch_present"),
  env_d = list(data = env_d, fixed = "site_in")
)

predictors <- c("cond", "depth", "oxy", "ph", "temp")
model_list <- list()

for (dname in names(datasets)) {
  dat          <- datasets[[dname]]$data
  fixed_effect <- datasets[[dname]]$fixed

  model_list[[dname]] <- list()

  for (resp in predictors) {
    form <- as.formula(paste(resp, "~", fixed_effect, "+ (1 | site_name) + (1 | sampling_round)"))
    model_list[[dname]][[resp]] <- lmer(form, data = dat, REML = TRUE)
  }
}


# --- 7.1 Model summary tables and coefficient plots --------------------------

tab_models_list  <- list()
plot_models_list <- list()

for (dname in names(model_list)) {
  models <- model_list[[dname]]

  tab_models_list[[dname]] <- tab_model(
    models[["cond"]], models[["depth"]], models[["oxy"]],
    models[["ph"]],   models[["temp"]],
    show.ngroups = FALSE,
    show.icc     = FALSE,
    show.zeroinf = FALSE,
    show.obs     = FALSE,
    show.se      = TRUE,
    collapse.ci  = TRUE,
    dv.labels    = c("Cond", "Depth", "Oxy", "pH", "Temp")
  )

  plot_models_list[[dname]] <- plot_models(
    models[["cond"]], models[["depth"]], models[["oxy"]],
    models[["ph"]],   models[["temp"]],
    show.p      = TRUE,
    show.values = TRUE,
    value.size  = 5
  ) + theme_pubr()

  print(plot_models_list[[dname]])
}

# Add descriptive titles to coefficient plots
plot_models_list[["env_m"]] <- plot_models_list[["env_m"]] + ggtitle("Merged dataset (env_m)")
plot_models_list[["env_r"]] <- plot_models_list[["env_r"]] + ggtitle("Center only (env_r)")
plot_models_list[["env_d"]] <- plot_models_list[["env_d"]] + ggtitle("Fields with ditch (env_d)")

# Combined coefficient plot across all three datasets
(plot_models_list[["env_m"]] /
    plot_models_list[["env_r"]] /
    plot_models_list[["env_d"]]) +
  plot_layout(guides = "collect") &
  theme(legend.box = "vertical", legend.margin = margin(), legend.position = "top")


# --- 7.2 ANOVA tables with Benjamini–Hochberg correction ---------------------

pretty_names <- c(
  cond  = "Conductivity",
  oxy   = "Dissolved oxygen",
  ph    = "pH",
  temp  = "Water temperature",
  depth = "Water depth"
)

extract_anova <- function(model, variable_name, dataset_name) {
  an    <- anova(model)
  raw_p <- an[["Pr(>F)"]]
  tibble(
    Dataset   = dataset_name,
    Variable  = pretty_names[[variable_name]],
    `Sum Sq`  = an[["Sum Sq"]],
    `num.df`  = an[["NumDF"]],
    `den.df`  = an[["DenDF"]],
    statistic = an[["F value"]],
    p         = raw_p
  )
}

# Build combined ANOVA table across all datasets and variables
anova_raw <- map_dfr(
  names(model_list),
  function(dataset) {
    map_dfr(
      names(model_list[[dataset]]),
      function(var_name) {
        extract_anova(model_list[[dataset]][[var_name]], var_name, dataset)
      }
    )
  }
)

# Apply BH correction within each dataset
anova_results <- anova_raw %>%
  group_by(Dataset) %>%
  mutate(
    p_adj               = p.adjust(p, method = "BH"),
    `p (formatted)`     = format.pval(p,     digits = 3, eps = 0.001),
    `p_adj (formatted)` = format.pval(p_adj, digits = 3, eps = 0.001)
  ) %>%
  ungroup()

anova_results

