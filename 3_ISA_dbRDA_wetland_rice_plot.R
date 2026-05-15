# =============================================================================
# Script:  3_ISA_wetland_rice_plot.R
# Purpose: Indicator species analysis (ISA) comparing rice paddy fields and
#          reference wetlands, followed by distance-based redundancy analysis
#          (dbRDA) to identify environmental drivers of macroinvertebrate
#          community composition. Produces the ISA forest plot and the dbRDA
#          ordination panels for the manuscript.
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
library(patchwork)
library(readxl)
library(data.table)
library(indicspecies)
library(vegan)

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

subsample_area_m2 <- 0.25 * 0.25

# Pivot to long format and calculate density (ind/m²)
comm_data_long <- comm_data %>%
  pivot_longer(
    cols      = Coenangrionidae:Sialis,
    names_to  = "species",
    values_to = "count"
  ) %>%
  mutate(density = (count / total_subamples) / subsample_area_m2)

# Taxonomic group lookup for each species
species_to_group <- c(
  "Coenangrionidae"          = "Odonata",
  "Coenagrion_sp"            = "Odonata",
  "Enallagma_cyathigerum"    = "Odonata",
  "Ischnura_sp"              = "Odonata",
  "Aeshna_sp"                = "Odonata",
  "Anax_sp"                  = "Odonata",
  "Anax_epiphiger"           = "Odonata",
  "Anax_imperator"           = "Odonata",
  "Anax_parthenope"          = "Odonata",
  "Libellulinae"             = "Odonata",
  "Libellula_depressa"       = "Odonata",
  "Libellula_quadrimaculata" = "Odonata",
  "Orthetrum_sp"             = "Odonata",
  "Orthetrum_albistylum"     = "Odonata",
  "Orthetrum_brunneum"       = "Odonata",
  "Orthetrum_cancellatum"    = "Odonata",
  "Crocothemis_erythraea"    = "Odonata",
  "Sympetrum_sp"             = "Odonata",
  "Sympetrum_depressiusculum"= "Odonata",
  "Sympetrum_fonscolombii"   = "Odonata",
  "Sympetrum_striolatum"     = "Odonata",
  "Anisus"                   = "Gastropoda",
  "Bithynia"                 = "Gastropoda",
  "Gyraulus"                 = "Gastropoda",
  "Stagnicola"               = "Gastropoda",
  "Physa"                    = "Gastropoda",
  "Planorbis"                = "Gastropoda",
  "Radix"                    = "Gastropoda",
  "Lithoglyphus_naticoides"  = "Gastropoda",
  "Pisidium"                 = "Bivalvia",
  "Sphaerium"                = "Bivalvia",
  "Gammarus"                 = "Amphipod",
  "Asellus_aquaticus"        = "Isopoda",
  "Cloeon"                   = "Ephemeroptera",
  "Caenis"                   = "Ephemeroptera",
  "other_Adephaga"           = "Coleoptera",
  "Laccophilus"              = "Coleoptera",
  "Hydroporinae"             = "Coleoptera",
  "other_dytiscidae"         = "Coleoptera",
  "other_Polyphaga"          = "Coleoptera",
  "Hydrophilidae"            = "Coleoptera",
  "Tipuloidea"               = "Diptera",
  "Psychodidae"              = "Diptera",
  "Chaoboridae"              = "Diptera",
  "Dixidae"                  = "Diptera",
  "Anophelinae"              = "Diptera",
  "Culicinae"                = "Diptera",
  "Ceratopogonidae"          = "Diptera",
  "Orthocladiinae"           = "Diptera",
  "Tanypodinae"              = "Diptera",
  "Chironomini"              = "Diptera",
  "Tanytarsini"              = "Diptera",
  "Brachyceres"              = "Diptera",
  "Gerromorpha"              = "Hemiptera",
  "Corixidae"                = "Hemiptera",
  "Micronecta"               = "Hemiptera",
  "Notonectoidea"            = "Hemiptera",
  "Plea_minutissima"         = "Hemiptera",
  "Trichoptera"              = "Trichoptera",
  "Sialis"                   = "Megaloptera"
)

comm_data_long <- comm_data_long %>%
  mutate(sp_group = species_to_group[species])

# Aggregate density to site × habitat × year × sampling round level
sp_data_m1 <- comm_data_long %>%
  group_by(site_name, env, year, sampling_round, species) %>%
  summarise(
    count                   = mean(density),
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

# Environmental metadata with log(x + 1) transformation and z-score scaling
env_m <- sp_data_m1[, 1:16]
env_m[, c(8, 10:12, 14)] <- log1p(env_m[, c(8, 10:12, 14)])
env_m[, c(8, 10:14)]     <- scale(env_m[, c(8, 10:14)], center = TRUE)

env_m$env            <- as.factor(env_m$env)
env_m$year           <- as.factor(env_m$year)
env_m$sampling_round <- as.factor(env_m$sampling_round)
env_m$ch_NS          <- as.factor(env_m$ch_NS)

env_only <- env_m[, c(8, 10:14)]

# Species matrix
sp_data_m <- sp_data_m1[, 29:ncol(sp_data_m1)]


# =============================================================================
# 4. Indicator species analysis (ISA)
# =============================================================================

# --- 4.1 Association strength (point-biserial r) with bootstrap CIs ----------

set.seed(1)
res <- strassoc(sp_data_m, sp_data_m1$env, func = "r.g", nboot.ci = 1000)

# --- 4.2 Permutation-based p-values (stratified within site) -----------------

set.seed(1)
res_sig <- signassoc(
  sp_data_m,
  cluster = sp_data_m1$env,
  control = how(
    blocks = sp_data_m1$site_name,
    within = Within(type = "free"),
    nperm  = 9999
  )
)

# Convert to data.table and add corrections
res_sig_df         <- as.data.table(res_sig, keep.rownames = TRUE)
res_sig_df$species <- res_sig_df$rn
res_sig_df$pholm   <- p.adjust(res_sig_df$psidak, method = "holm")
res_sig_df$pfdr    <- p.adjust(res_sig_df$psidak, method = "BH")

# Significant species under different correction schemes
res_sig_df[pfdr  <= 0.05, ]
res_sig_df[pholm <= 0.05, ]
res_sig_df[psidak <= 0.05, ]

# Assign habitat direction
res_sig_df[, habitat := fifelse(best == 1, "rice",
                                fifelse(best == 2, "wetland", NA_character_))]


# --- 4.3 Build plot data frame -----------------------------------------------

rice_col    <- "tomato2"
wetland_col <- "#0072b2"

stat  <- as.data.frame(res$stat)
lower <- as.data.frame(res$lowerCI)
upper <- as.data.frame(res$upperCI)

assoc              <- stat$rice
assoc[is.na(assoc)] <- 0

df_plot <- data.frame(
  species  = rownames(stat),
  assoc    = assoc,
  lowerCI  = lower$rice,
  upperCI  = upper$rice
) %>%
  left_join(res_sig_df, by = c("species" = "rn")) %>%
  mutate(
    assoc_dir = ifelse(assoc > 0, "Rice", "Wetland"),
    sig = case_when(
      pfdr <= 0.001 ~ "***",
      pfdr <= 0.01  ~ "**",
      pfdr <= 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    group = species_to_group[species]
  )

# Taxonomic group order
group_order  <- c("Gastropoda", "Bivalvia", "Amphipod", "Isopoda",
                  "Ephemeroptera", "Odonata", "Megaloptera",
                  "Trichoptera", "Coleoptera", "Diptera", "Hemiptera")
df_plot$group <- factor(df_plot$group, levels = group_order)

# Species order within groups (sorted by association strength)
species_order <- df_plot %>%
  arrange(group, assoc) %>%
  pull(species)

df_plot <- df_plot %>%
  mutate(
    species = factor(species, levels = species_order),
    y_pos   = as.numeric(species),
    label_color = case_when(
      pfdr > 0.05          ~ "black",
      assoc_dir == "Rice"    ~ rice_col,
      assoc_dir == "Wetland" ~ wetland_col,
      TRUE                   ~ "black"
    )
  )

# Group annotation bar positions
group_lines <- df_plot %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(
    y_min = min(y_pos),
    y_max = max(y_pos),
    y_mid = (y_min + y_max) / 2,
    .groups = "drop"
  )

label_x <- min(df_plot$assoc, na.rm = TRUE) -
  0.3 * (max(df_plot$assoc, na.rm = TRUE) - min(df_plot$assoc, na.rm = TRUE))
xmax <- max(df_plot$assoc, na.rm = TRUE)
xmin <- min(df_plot$assoc, na.rm = TRUE)

# Long-format data frame with FDR significance and color assignment
df_long <- df_plot %>%
  dplyr::select(species, assoc, lowerCI, upperCI, assoc_dir, pfdr, habitat, group, y_pos) %>%
  mutate(
    signif = case_when(
      pfdr <= 0.001 ~ "***",
      pfdr <= 0.01  ~ "**",
      pfdr <= 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    assoc_pos = assoc,
    color_env = case_when(
      signif == "ns"         ~ "black",
      habitat == "rice"      ~ rice_col,
      habitat == "wetland"   ~ wetland_col,
      TRUE                   ~ "black"
    )
  )


# --- 4.4 ISA forest plot -----------------------------------------------------

p_indicator3_grouped <- ggplot(df_plot, aes(x = assoc, y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI),
                 height = 0.2, color = "gray50", size = 0.8) +
  # Significance symbol per species
  geom_point(data = df_long,
             aes(x = assoc_pos, y = species, shape = signif, color = color_env),
             size = 3.5, show.legend = FALSE) +
  # Species labels (wetland side)
  geom_text(data = df_plot %>% filter(assoc_dir == "Wetland"),
            aes(x = lowerCI, y = species,
                label = gsub("_", " ", species), color = label_color),
            hjust = 1.1, vjust = 0.3, size = 4,
            fontface = "bold.italic", show.legend = FALSE) +
  # Species labels (rice side)
  geom_text(data = df_plot %>% filter(assoc_dir == "Rice"),
            aes(x = upperCI, y = species,
                label = gsub("_", " ", species), color = label_color),
            hjust = -0.1, vjust = 0.3, size = 4,
            fontface = "bold.italic", show.legend = FALSE) +
  # Taxonomic group annotation bars
  geom_segment(data = group_lines,
               aes(x = label_x - 2, xend = label_x - 2, y = y_min, yend = y_max),
               color = "gray30", size = 1, inherit.aes = FALSE) +
  geom_text(data = group_lines,
            aes(x = label_x - 2.05, y = y_mid, label = group),
            hjust = 1, fontface = "bold", color = "gray30",
            size = 4, inherit.aes = FALSE) +
  scale_shape_manual(values = c("ns" = 39, "*" = 15, "**" = 17, "***" = 16)) +
  scale_color_identity() +
  # Direction arrows
  annotate("segment", x = 0, xend = xmax * 0.9,
           y = length(levels(df_plot$species)) + 4,
           yend = length(levels(df_plot$species)) + 4,
           arrow = arrow(length = unit(0.5, "cm")),
           color = rice_col, size = 3) +
  annotate("text", x = xmax * 1.0,
           y = length(levels(df_plot$species)) + 4.5,
           label = "Rice", color = rice_col,
           hjust = 0, vjust = 0, size = 4.5) +
  annotate("segment", x = 0, xend = xmin * 0.9,
           y = length(levels(df_plot$species)) + 4,
           yend = length(levels(df_plot$species)) + 4,
           arrow = arrow(length = unit(0.5, "cm")),
           color = wetland_col, size = 3) +
  annotate("text", x = xmin * 1.0,
           y = length(levels(df_plot$species)) + 4.5,
           label = "Wetland", color = wetland_col,
           hjust = 1, vjust = 0, size = 4.5) +
  labs(x = "Point-biserial correlation coefficient", y = NULL,
       shape = "Significance level") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    panel.grid.major.y      = element_blank(),
    axis.text.y             = element_blank(),
    axis.ticks.y            = element_blank(),
    legend.position         = "bottom",
    legend.box              = "vertical",
    plot.margin             = margin(10, 100, 10, 100, unit = "pt"),
    legend.background       = element_rect(fill = "white", colour = "white")
  )

p_indicator3_grouped


# =============================================================================
# 5. Distance-based Redundancy Analysis (dbRDA)
# =============================================================================

# Contrasts set for balanced Type III ANOVA tests
options(contrasts = c("contr.sum", "contr.poly"))

# Square-root transformed species matrix and Bray-Curtis dissimilarity matrix
hell_m  <- sqrt(sp_data_m1[, 29:ncol(sp_data_m1)])
brayS_m <- sqrt(vegdist(hell_m, method = "bray"))


# --- 5.1 Homogeneity of multivariate dispersion ------------------------------

beta_disp_m <- betadisper(brayS_m, group = env_m$env)
plot(beta_disp_m)

permutest(beta_disp_m, permutations = how(
  blocks = env_m$site_name,
  within = Within(type = "free"),
  nperm  = 999
))


# --- 5.2 Forward selection of environmental predictors -----------------------

# Full model including habitat type and abiotic variables; site identity, year,
# and sampling round included as conditioning variables
dbRDA_env <- dbrda(
  brayS_m ~ env + area_ha + cond + ph + oxy + temp + depth +
    Condition(site_name, year, sampling_round),
  data = env_m
)

# Null model (conditioning variables only)
dbRDA_m0 <- dbrda(
  brayS_m ~ 1 + Condition(site_name, year, sampling_round),
  data = env_m
)

adjR2.dbRDA_env <- RsquareAdj(dbRDA_env)$adj.r.squared

set.seed(1)
ordi_p_env <- ordiR2step(
  dbRDA_m0,
  scope        = formula(dbRDA_env),
  R2scope      = adjR2.dbRDA_env,
  direction    = "forward",
  permutations = 999
)
ordi_p_env$anova

# Apply BH correction to forward selection p-values
ordi_p_env_adj <- ordi_p_env
ordi_p_env_adj$anova$`Pr(>F)` <- p.adjust(
  ordi_p_env$anova$`Pr(>F)`,
  method = "BH",
  n      = ncol(env_m)
)
ordi_p_env_adj$anova


# --- 5.3 Permutation tests on the selected model -----------------------------

# Permutation design: free within site, blocked by site identity
perm_block <- how(
  blocks = env_m$site_name,
  within = Within(type = "free"),
  nperm  = 999
)

set.seed(1)
glob_env.aov <- anova.cca(ordi_p_env_adj, permutations = perm_block)

set.seed(1)
terms_env.aov <- anova.cca(ordi_p_env_adj, permutations = perm_block, by = "terms")
terms_env.aov$`Pr(>F)` <- p.adjust(terms_env.aov$`Pr(>F)`, method = "BH")

set.seed(1)
axis_env.aov <- anova.cca(ordi_p_env_adj, permutations = perm_block, by = "axis")
axis_env.aov$`Pr(>F)` <- p.adjust(axis_env.aov$`Pr(>F)`, method = "BH")

summary(ordi_p_env_adj)
RsquareAdj(ordi_p_env_adj)


# --- 5.4 Extract ordination scores -------------------------------------------

# Add species scores from the Hellinger-transformed matrix
sppscores(ordi_p_env_adj) <- hell_m

scores_e <- scores(ordi_p_env_adj, choices = 1:3)
sp_e  <- as.data.frame(scores_e[["species"]][, 1:3])
var_e <- as.data.frame(scores_e[["biplot"]][, 1:3])
st_e  <- as.data.frame(scores_e[["sites"]][, 1:3])

# Retain only ISA-significant species for ordination overlay
sp_e$species <- rownames(sp_e)
sig_sp       <- df_plot %>% filter(pfdr <= 0.05)

sp_e <- sp_e %>%
  filter(species %in% sig_sp$species) %>%
  left_join(sig_sp[, c("species", "habitat")], by = "species")

# Add metadata to site scores
st_e          <- cbind(st_e, env_m[, 1:ncol(env_m)])
st_e$habitat  <- st_e$env
st_e$sampling_period <- factor(paste0(st_e$year, "_", st_e$sampling_round))

# Convex hulls per habitat
r_e <- st_e[st_e$env == "rice", ][chull(st_e[st_e$env == "rice",    c("dbRDA1", "dbRDA2")]), ]
w_e <- st_e[st_e$env == "wetland", ][chull(st_e[st_e$env == "wetland", c("dbRDA1", "dbRDA2")]), ]
hull.data_e         <- rbind(r_e, w_e)
hull.data_e$habitat <- hull.data_e$env


# --- 5.5 Axis labels and display settings ------------------------------------

sampling_labels <- c(
  "2022_1" = "July '22",
  "2022_2" = "August '22",
  "2023_1" = "July '23",
  "2023_2" = "August '23"
)

var_labels <- c(
  "env1"    = "Habitat type",
  "habitat" = "Habitat type",
  "cond"    = "Conductivity",
  "depth"   = "Depth",
  "oxy"     = "Dissolved \noxygen",
  "area_ha" = "Water surface \narea"
)
rownames(var_e) <- var_labels[rownames(var_e)]

# Axis labels: % explained variance of constrained axes (fitted) and total inertia
pct_fitted <- 100 * ordi_p_env_adj$CCA$eig / ordi_p_env_adj$CCA$tot.chi
pct_total  <- 100 * summary(ordi_p_env_adj)$cont$importance[2, ]

x_lab <- paste0("dbRDA1 (", round(pct_fitted[1], 1), "% fitted; ",
                round(pct_total[1], 1), "% total)")
y_lab <- paste0("dbRDA2 (", round(pct_fitted[2], 1), "% fitted; ",
                round(pct_total[2], 1), "% total)")


# --- 5.6 dbRDA plots ---------------------------------------------------------

# Panel A: Site scores coloured by habitat with environmental vectors
p_e <- ggplot() +
  geom_polygon(data = hull.data_e,
               aes(x = dbRDA1, y = dbRDA2, fill = habitat, colour = habitat),
               alpha = 0.1, show.legend = FALSE) +
  geom_point(data = st_e,
             aes(x = dbRDA1, y = dbRDA2, shape = sampling_period, color = habitat),
             size = 2.5, alpha = 0.5) +
  geom_segment(data = var_e,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "gray20", linewidth = 1) +
  geom_text_repel(data = var_e,
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(var_e)),
                  size = 5, box.padding = 0.8, color = "gray20",
                  force = 2, force_pull = 5, fontface = "bold") +
  scale_shape_manual(values = c(16, 1, 17, 2, 15, 18),
                     labels = sampling_labels,
                     name   = "Sampling month & year") +
  scale_color_manual(values = c(rice = rice_col, wetland = wetland_col),
                     name   = "Habitat") +
  scale_fill_manual(values  = c(rice = rice_col, wetland = wetland_col)) +
  labs(x = x_lab, y = y_lab, tag = "A") +
  guides(
    shape  = guide_legend(nrow = 2, order = 1),
    colour = guide_legend(nrow = 2, order = 1)
  ) +
  theme_pubr() +
  theme(legend.position = "top")

# Panel B: Significant species scores overlaid on habitat convex hulls
sp_e$species <- gsub("_", " ", sp_e$species)

p2e <- ggplot() +
  geom_polygon(data = hull.data_e,
               aes(x = dbRDA1, y = dbRDA2, fill = habitat, colour = habitat),
               alpha = 0.1, show.legend = FALSE) +
  geom_point(data = sp_e,
             aes(x = dbRDA1, y = dbRDA2),
             fill = "black", color = "grey", shape = 22, size = 2) +
  geom_text_repel(data = sp_e,
                  aes(x = dbRDA1, y = dbRDA2, label = species, color = habitat),
                  box.padding = 0.5, max.overlaps = Inf,
                  size = 4, force = 2, fontface = "bold.italic",
                  show.legend = FALSE) +
  scale_color_manual(values = c(rice = rice_col, wetland = wetland_col),
                     name   = "Habitat") +
  scale_fill_manual(values  = c(rice = rice_col, wetland = wetland_col)) +
  labs(x = x_lab, y = y_lab, tag = "B") +
  theme_pubr()

# Combined panels A + B
p_e / p2e +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.title.position = "top")


# --- 5.7 Combined triplot: sites + significant species + environment ----------

# Offset label positions slightly beyond arrow tips
var_e$lab_x <- var_e$dbRDA1 + 1.5
var_e$lab_y <- var_e$dbRDA2 + 1.5

p_e2 <- ggplot() +
  geom_polygon(data = hull.data_e,
               aes(x = dbRDA1, y = dbRDA2, fill = habitat, colour = habitat),
               alpha = 0.1, show.legend = FALSE) +
  geom_point(data = st_e,
             aes(x = dbRDA1, y = dbRDA2, shape = sampling_period, color = habitat),
             size = 2.5, alpha = 0.5) +
  geom_point(data = sp_e,
             aes(x = dbRDA1, y = dbRDA2, fill = habitat),
             shape = 22, size = 3, color = "black", show.legend = FALSE) +
  geom_label_repel(data = sp_e,
                   aes(x     = dbRDA1, y = dbRDA2,
                       label = species, color = habitat,
                       fill  = after_scale(scales::alpha(fill, 0.6))),
                   box.padding = 0.5, max.overlaps = Inf,
                   size = 4, force = 5, fontface = "bold.italic",
                   show.legend = FALSE, segment.size = 1) +
  geom_segment(data = var_e,
               aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "gray20", linewidth = 1.5) +
  geom_text_repel(data = var_e,
                  aes(x = dbRDA1, y = dbRDA2, label = rownames(var_e)),
                  size = 5, color = "gray20", fontface = "bold",
                  box.padding = 0.3, point.padding = 0.2,
                  segment.color = NA, max.overlaps = Inf) +
  scale_shape_manual(values = c(16, 1, 17, 2, 15, 18),
                     labels = sampling_labels,
                     name   = "Sampling month & year") +
  scale_color_manual(values = c(rice = rice_col, wetland = wetland_col),
                     name   = "Habitat") +
  scale_fill_manual(values  = c(rice = rice_col, wetland = wetland_col)) +
  labs(x = x_lab, y = y_lab) +
  guides(
    shape  = guide_legend(nrow = 2, order = 1),
    colour = guide_legend(nrow = 2, order = 1)
  ) +
  theme_pubr() +
  theme(legend.position = "top")

p_e2
