# =============================================================================
# Script:  30.07.2025_ISA_rice_field_comparisons.R
# Purpose: Perform Indicator Species Analysis (ISA) across three datasets:
#          (1) rice field center vs. ditch samples (only sites with ditches),
#          (2) merged dataset comparing fields with vs. without ditches,
#          (3) rice field center-only dataset comparing fields with vs. without ditches.
#          Includes effect size estimation, permutation-based significance testing,
#          multiple-testing correction, and visualization of indicator species.
# Author:  Thea Bulas
# Manuscript title: Rice paddies promote diverse and distinct aquatic
#                   invertebrate communities in agroecosystems.
# Date:    15.05.2026
# Data:    2022_2023_macroinvRADICAL_abio_env_20250407.xlsx
# =============================================================================



#library(vegan3d)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggpp)
library(ggrepel)
library(patchwork)
library(readxl)
library(PerformanceAnalytics)
library(corrplot)
library(indicspecies)
library(ggh4x)
library(lme4)
library(lmerTest)
library(data.table)




setwd("~/RStudio/agroscope/chapter_2/Manuscript_scripts")


# -----------------------------------------------------------------------------
# Data import
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




#add info about ditch presence:
# Define site names with ditch presence
ditch_yes <- c("brugg", "detligen", "jonen", "muehlau", "untersiggenthal_ost", "vionnaz")
ditch_no <- c("la_sauge", "stetten", "untersiggenthal_west", "witzwil", "wuerenlingen")

# Add 'ditch_present' column to comm_data
comm_data1 <- comm_data %>%
  mutate(
    ditch_present = case_when(
      site_name %in% ditch_yes ~ "yes",
      site_name %in% ditch_no ~ "no",
      TRUE ~ NA_character_
    ),
    .after = site_name
  )

comm_data_long<- comm_data1 %>% filter(env == "rice") %>% pivot_longer(cols=c("Coenangrionidae":"Sialis"), names_to='species', values_to='count')


#### FULL sp count matrix (rice CENTER and DITCH kept separate)
subsample_area_m2 <- (0.25 * 0.25)

sp_data_f <- comm_data_long %>%
  group_by(site_name, site_in, year, sampling_round, species) %>%
  summarise(count = mean(count/total_subamples)/subsample_area_m2, #to get ind/m2
            total_subsamples=sum(total_subamples),
            ditch_present=first(ditch_present),
            canton=first(canton),
            ch_NS=first(ch_NS),
            area_ha=sum(area_ha),
            distance_rice_wetland_m = mean(distance_rice_wetland_m),
            #abio measurments
            cond=mean(cond),
            depth=mean(depth),
            oxy=mean(oxy),
            ph=mean(ph),
            temp=mean(temp),
            NO=mean(NO),
            PO=mean(PO),
            
            
            #fertiliser:
            N_fertiliser_name=first(N_fertiliser_name),
            fertiliser_type=first(fertiliser_type),
            total_N_per_ha=mean(total_N_per_ha),
            fertiliser_latest_date = min(fertiliser_latest_date),
            #soils:
            soil_prep_type=first(soil_prep_type),
            soil_prep_depth_cm=mean(soil_prep_depth_cm),
            soil_prep_date=min(soil_prep_date),
            #rice plant in and out
            rice_transplant=min(rice_transplant),
            rice_harvest=min(rice_harvest),
            #water managment
            water_overwinter=min(water_overwinter),
            water_start=min(water_start),
            water_stop=min(water_stop)) %>%
  ungroup() %>%
  pivot_wider(names_from="species", values_from="count")

meta_cols_f <- c("site_name", "site_in", "year", "sampling_round", 
                 "total_subsamples", "ditch_present", "canton", "ch_NS", 
                 "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy", 
                 "ph", "temp", "NO", "PO", "N_fertiliser_name", 
                 "fertiliser_type", "total_N_per_ha", "fertiliser_latest_date", 
                 "soil_prep_type", "soil_prep_depth_cm", "soil_prep_date", 
                 "rice_transplant", "rice_harvest", "water_overwinter", 
                 "water_start", "water_stop")


# Remove species columns with only 0 or NA values
#if issue with not working argument, try to runn dply::select. SOme functions may be not working due to other packages
sp_data_f_cleaned <- sp_data_f %>%
  dplyr::select(all_of(meta_cols_f), where(~ any(. != 0, na.rm = TRUE)))


###
### Dataset for only fields with ditch:
sp_data_ditch <- comm_data_long %>%
  filter(ditch_present == "yes") %>%
  group_by(site_name, site_in, year, sampling_round, species) %>%
  summarise(count = mean(count/total_subamples)/subsample_area_m2, #to get ind/m2
            total_subsamples=sum(total_subamples),
            ditch_present=first(ditch_present),
            canton=first(canton),
            ch_NS=first(ch_NS),
            area_ha=sum(area_ha),
            distance_rice_wetland_m = mean(distance_rice_wetland_m),
            #abio measurments
            cond=mean(cond),
            depth=mean(depth),
            oxy=mean(oxy),
            ph=mean(ph),
            temp=mean(temp),
            NO=mean(NO),
            PO=mean(PO),
            
            
            #fertiliser:
            N_fertiliser_name=first(N_fertiliser_name),
            fertiliser_type=first(fertiliser_type),
            total_N_per_ha=mean(total_N_per_ha),
            fertiliser_latest_date = min(fertiliser_latest_date),
            #soils:
            soil_prep_type=first(soil_prep_type),
            soil_prep_depth_cm=mean(soil_prep_depth_cm),
            soil_prep_date=min(soil_prep_date),
            #rice plant in and out
            rice_transplant=min(rice_transplant),
            rice_harvest=min(rice_harvest),
            #water managment
            water_overwinter=min(water_overwinter),
            water_start=min(water_start),
            water_stop=min(water_stop)) %>%
  ungroup() %>%
  pivot_wider(names_from="species", values_from="count")

meta_cols <- c("site_name", "site_in", "year", "sampling_round", 
               "total_subsamples", "ditch_present", "canton", "ch_NS", 
               "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy", 
               "ph", "temp", "NO", "PO", "N_fertiliser_name", 
               "fertiliser_type", "total_N_per_ha", "fertiliser_latest_date", 
               "soil_prep_type", "soil_prep_depth_cm", "soil_prep_date", 
               "rice_transplant", "rice_harvest", "water_overwinter", 
               "water_start", "water_stop")

# Remove species columns with only 0 or NA values
sp_data_ditch_cleaned <- sp_data_ditch %>%
  dplyr::select(all_of(meta_cols), where(~ any(. != 0, na.rm = TRUE)))

sp_data_d <- sp_data_ditch_cleaned[,(30:ncol(sp_data_ditch_cleaned))]
env_d <- sp_data_ditch_cleaned[,(1:29)]
env_d[,c(9,11:13,15)] <- log1p(env_d[,c(9,11:13,15)])
env_d[,c(9,11:15)] <- scale(env_d[,c(9,11:15)])

env_only_d <- env_d[,c(9,11:15)]



#
#Create named vector for species-group mapping
species_to_group <- c(
  "Coenangrionidae" = "Odonata",
  "Coenagrion_sp" = "Odonata",
  "Enallagma_cyathigerum" = "Odonata",
  "Ischnura_sp" = "Odonata",
  "Aeshna_sp" = "Odonata",
  "Anax_sp" = "Odonata",
  "Anax_epiphiger" = "Odonata",
  "Anax_imperator" = "Odonata",
  "Anax_parthenope" = "Odonata",
  "Libellulinae" = "Odonata",
  "Libellula_depressa" = "Odonata",
  "Libellula_quadrimaculata" = "Odonata",
  "Orthetrum_sp" = "Odonata",
  "Orthetrum_albistylum" = "Odonata",
  "Orthetrum_brunneum" = "Odonata",
  "Orthetrum_cancellatum" = "Odonata",
  "Crocothemis_erythraea" = "Odonata",
  "Sympetrum_sp" = "Odonata",
  "Sympetrum_depressiusculum" = "Odonata",
  "Sympetrum_fonscolombii" = "Odonata",
  "Sympetrum_striolatum" = "Odonata",
  "Anisus" = "Gastropoda",
  "Bithynia" = "Gastropoda",
  "Gyraulus" = "Gastropoda",
  "Stagnicola" = "Gastropoda",
  "Physa" = "Gastropoda",
  "Planorbis" = "Gastropoda",
  "Radix" = "Gastropoda",
  "Lithoglyphus_naticoides" = "Gastropoda",
  "Pisidium" = "Bivalvia",
  "Sphaerium" = "Bivalvia",
  "Gammarus" = "Amphipod",
  "Asellus_aquaticus" = "Isopoda",
  "Cloeon" = "Ephemeroptera",
  "Caenis" = "Ephemeroptera",
  "other_Adephaga" = "Coleoptera",
  "Laccophilus" = "Coleoptera",
  "Hydroporinae" = "Coleoptera",
  "other_dytiscidae" = "Coleoptera",
  "other_Polyphaga" = "Coleoptera",
  "Hydrophilidae" = "Coleoptera",
  "Tipuloidea" = "Diptera",
  "Psychodidae" = "Diptera",
  "Chaoboridae" = "Diptera",
  "Dixidae" = "Diptera",
  "Anophelinae" = "Diptera",
  "Culicinae" = "Diptera",
  "Ceratopogonidae" = "Diptera",
  "Orthocladiinae" = "Diptera",
  "Tanypodinae" = "Diptera",
  "Chironomini" = "Diptera",
  "Tanytarsini" = "Diptera",
  "Brachyceres" = "Diptera",
  "Gerromorpha" = "Hemiptera",
  "Corixidae" = "Hemiptera",
  "Micronecta" = "Hemiptera",
  "Notonectoidea" = "Hemiptera",
  "Plea_minutissima" = "Hemiptera",
  "Trichoptera" = "Trichoptera",
  "Sialis" = "Megaloptera"
)



#

#calculate effect size and CI
set.seed(1)
res <- strassoc(sp_data_d, sp_data_ditch_cleaned$site_in, func="r.g", nboot.ci = 1000)
res

# permutation-based p-values (with stratification)
set.seed(1)
res_sig <- signassoc(sp_data_d, cluster=sp_data_ditch_cleaned$site_in, mode = 1, control = how(
  blocks = sp_data_ditch_cleaned$site_name,  # Only restrict permutations to within-site
  within = Within(type = "free"),
  nperm = 9999))

res_sig_df <- as.data.table(res_sig, keep.rownames=TRUE)
res_sig_df$species <- res_sig_df$rn

res_sig_df$pfdr <- p.adjust(res_sig_df$psidak, method = "BH")

res_sig_df[pfdr<=0.05, ]


#Add habitat information automatically
res_sig_df[, habitat := fifelse(best == 1, "center",
                                fifelse(best == 2, "side", NA_character_))]

# Extract results
stat <- as.data.frame(res$stat)
lower <- as.data.frame(res$lowerCI)
upper <- as.data.frame(res$upperCI)

# Compute the rice-wetland difference
assoc <- stat$center  # This is positive toward rice
assoc[is.na(assoc)] <- 0

# Combine into one data frame
df_plot_d <- data.frame(
  species = rownames(stat),
  assoc = assoc,
  lowerCI = lower$center,
  upperCI = upper$center
)

# Merge with significance info
df_plot_d <- left_join(df_plot_d, res_sig_df, by = c("species" = "rn"))

# Add direction (rice vs wetland)
df_plot_d <- df_plot_d %>%
  mutate(
    assoc_dir = ifelse(assoc > 0, "Field", "Ditch"),
    sig_fdr = case_when(
      pfdr <= 0.001 ~ "***",
      pfdr <= 0.01 ~ "**",
      pfdr <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

# Sort by association strength
#df_plot$species <- factor(df_plot$species, levels = df_plot$species[order(df_plot$assoc)])



library(wesanderson)
wes_col <- wes_palette("FantasticFox1", type = "discrete")
center_col <- "orange3"
side_col <- "gold3"

# 1. Reshape into long format
df_long_d <- df_plot_d %>%
  dplyr::select(species, assoc, lowerCI, upperCI, assoc_dir, psidak, sig_fdr) %>%
  pivot_longer(cols = c(sig_fdr),
               names_to = "correction", values_to = "signif") %>%
  mutate(
    correction = case_when(
      correction == "psidak" ~ "Sidak",
      correction == "sig_fdr" ~ "FDR"
    ),
    # x-offset to avoid overlap: -0.015, 0, +0.015
    x_offset = case_when(
      correction == "Sidak" ~ -0.05,
      correction == "FDR"   ~ 0.05
    ),
    assoc_pos = assoc + x_offset
  )

df_long_d <- df_long_d %>%
  mutate(sp_group = species_to_group[species])

# Define taxonomic chronological order of groups
group_order <- c("Gastropoda", "Bivalvia", "Amphipod", "Isopoda",
                 "Ephemeroptera", "Odonata", "Megaloptera",
                 "Trichoptera", "Coleoptera", "Diptera", "Hemiptera")


# Reorder species within the plot by group first, then assoc
df_plot_d <- df_plot_d %>%
  mutate(group = factor(species_to_group[species], levels = group_order)) %>% # preserve order
  arrange(group, desc(assoc))

# Reset species as factor with the new order
df_plot_d$species <- factor(df_plot_d$species, levels = df_plot_d$species)
# Apply the same ordering to df_long_d
df_long_d$species <- factor(df_long_d$species, levels = levels(df_plot_d$species))

# Recompute y positions for vertical group lines
df_plot_d$y_pos <- as.numeric(factor(df_plot_d$species))

group_lines_d <- df_plot_d %>%
  group_by(group) %>%
  summarize(
    y_min = min(y_pos),
    y_max = max(y_pos),
    y_mid = (min(y_pos) + max(y_pos)) / 2
  )

# Set x position for group labels
label_x <- -1  # adjust left of the plot




##
p_indicator_ditch <- ggplot(df_plot_d, aes(x = assoc, y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.2,
                 color = "gray50", size = 0.8) +
  
  # Add three significance symbols per species
  geom_point(data = df_long_d,
             aes(x = assoc_pos, y = species,
                 shape = signif, color = correction),
             size = 3.5) +
  
  # Species labels
  # Field (positive assoc) = side_col
  # Species labels (separated by habitat side)
  geom_text(data = df_plot_d[df_plot_d$assoc_dir == "Field", ],
            aes(label = gsub("_", " ", species)),
            hjust = -0.1, vjust = 0.3,
            x = df_plot_d$upperCI[df_plot_d$assoc_dir == "Field"],
            size = 4, fontface = "bold.italic", color = side_col) +   # Field = side_col (gold3)
  
  geom_text(data = df_plot_d[df_plot_d$assoc_dir == "Ditch", ],
            aes(label = gsub("_", " ", species)),
            hjust = 1.1, vjust = 0.3,
            x = df_plot_d$lowerCI[df_plot_d$assoc_dir == "Ditch"],
            size = 4, fontface = "bold.italic", color = center_col) + # Ditch = center_col (orange3)
  
  
  # Group annotation lines
  geom_segment(data = group_lines_d,
               aes(x = label_x, xend = label_x, y = y_min, yend = y_max),
               color = "gray30", size = 1) +
  
  # Group labels (no exclamation marks)
  geom_text(data = group_lines_d,
            aes(x = label_x - 0.05, y = (y_min + y_max)/2, label = group),
            hjust = 1, fontface = "bold", color = "gray30", size = 4)+
  
  # Shapes & colors for significance
  scale_shape_manual(values = c("ns" = 39, "*" = 15, "**" = 17, "***" = 16)) +
  scale_color_manual(values = c("Sidak" = "gold3",
                                "Holm"  = "purple",
                                "FDR"   = "darkgreen")) +
  # Arrow + label for Field (right / positive)
  annotate("segment", x = 0, xend = max(df_plot_d$assoc) * 0.9, 
           y = length(unique(df_plot_d$species)) + 4, 
           yend = length(unique(df_plot_d$species)) + 4, 
           arrow = arrow(length = unit(0.5, "cm")), 
           color = side_col, size = 3) +
  annotate("text", x = max(df_plot_d$assoc) * 1, 
           y = length(unique(df_plot_d$species)) + 4.5, 
           label = "Field", color = side_col, hjust = 0, vjust = 0, size = 4.5) +
  
  # Arrow + label for Ditch (left / negative)
  annotate("segment", x = 0, xend = min(df_plot_d$assoc) * 0.9, 
           y = length(unique(df_plot_d$species)) + 4, 
           yend = length(unique(df_plot_d$species)) + 4, 
           arrow = arrow(length = unit(0.5, "cm")), 
           color = center_col, size = 3) +
  annotate("text", x = min(df_plot_d$assoc) * 1, 
           y = length(unique(df_plot_d$species)) + 4.5, 
           label = "Ditch", color = center_col, hjust = 1, vjust = 0, size = 4.5)+
  theme_minimal() +
  labs(
    x = "Point-Biserial Correlation Coefficient",
    y = NULL,
    color = "Correction method",
    shape = "Significance level"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 100, 10, 100, unit = "pt"),
    legend.background = element_rect(fill = "white", colour = "white")
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))

#FIgure 8
p_indicator_ditch



###### Ditch present absent overall between rice fields with:
## MErged dataset

#### MERGED sp count matrix
#avaraged across site_in and pond_# and calculate densities (ind/m2)

sp_data_m1 <- comm_data_long %>%
  group_by(site_name, year, sampling_round, species) %>% #we don't include "site_in" here to aggregate
  summarise(count = mean(count/total_subamples)/subsample_area_m2, #to get ind/m2
            total_subsamples=sum(total_subamples),
            ditch_present=first(ditch_present),
            canton=first(canton),
            ch_NS=first(ch_NS),
            area_ha=sum(area_ha),
            distance_rice_wetland_m = mean(distance_rice_wetland_m),
            #abio measurments
            cond=mean(cond),
            depth=mean(depth),
            oxy=mean(oxy),
            ph=mean(ph),
            temp=mean(temp),
            NO=mean(NO),
            PO=mean(PO),
            
            
            #fertiliser:
            N_fertiliser_name=first(N_fertiliser_name),
            fertiliser_type=first(fertiliser_type),
            total_N_per_ha=mean(total_N_per_ha),
            fertiliser_latest_date = min(fertiliser_latest_date),
            #soils:
            soil_prep_type=first(soil_prep_type),
            soil_prep_depth_cm=mean(soil_prep_depth_cm),
            soil_prep_date=min(soil_prep_date),
            #rice plant in and out
            rice_transplant=min(rice_transplant),
            rice_harvest=min(rice_harvest),
            #water managment
            water_overwinter=min(water_overwinter),
            water_start=min(water_start),
            water_stop=min(water_stop)) %>%
  ungroup() %>%
  pivot_wider(names_from="species", values_from="count")

meta_cols_m <- c("site_name", "year", "sampling_round", 
                 "total_subsamples", "ditch_present", "canton", "ch_NS", 
                 "area_ha", "distance_rice_wetland_m", "cond", "depth", "oxy", 
                 "ph", "temp", "NO", "PO", "N_fertiliser_name", 
                 "fertiliser_type", "total_N_per_ha", "fertiliser_latest_date", 
                 "soil_prep_type", "soil_prep_depth_cm", "soil_prep_date", 
                 "rice_transplant", "rice_harvest", "water_overwinter", 
                 "water_start", "water_stop")

# Remove species columns with only 0 or NA values
sp_data_m_cleaned <- sp_data_m1 %>%
  dplyr::select(all_of(meta_cols_m), where(~ any(. != 0, na.rm = TRUE)))

# MERGED (ditch+center) sp count matrix
sp_data_m <- sp_data_m_cleaned[,29:ncol(sp_data_m_cleaned)]
# Create env data frame and log+1 transform
env_m <- sp_data_m_cleaned[,1:28]
#Log+1 on cond, oxy, temp, depth, area_ha
env_m[,c(8,10:12,14)] <- log1p(env_m[,c(8,10:12,14)])
env_m[,c(8,10:14)] <- scale(env_m[,c(8,10:14)], center = T)


env_m$year <- as.factor(env_m$year)
env_m$sampling_round <- as.factor(env_m$sampling_round)
env_m$ch_NS<- as.factor(env_m$ch_NS)

#env only dataframe (no meta data)
env_m_only <- env_m[,c(8,10:14)]


##
#calculate effect size and CI
set.seed(1)
res_m <- strassoc(sp_data_m, sp_data_m_cleaned$ditch_present, func="r.g", nboot.ci = 1000)
res_m


# permutation-based p-values (with stratification)
set.seed(1)
res_sig_m <- signassoc(sp_data_m, cluster=sp_data_m_cleaned$ditch_present, mode = 1, control = how(
  blocks = sp_data_m_cleaned$site_name,  # Only restrict permutations to within-site
  within = Within(type = "free"),
  nperm = 9999))

res_sig_m_df <- as.data.table(res_sig_m, keep.rownames=TRUE)
res_sig_m_df$species <- res_sig_m_df$rn

#corrections
#res_sig_m_df$pholm <- p.adjust(res_sig_m_df$psidak, method = "holm")
res_sig_m_df$pfdr <- p.adjust(res_sig_m_df$psidak, method = "BH")

res_sig_m_df[psidak<=0.05, ]


#Add habitat information automatically
res_sig_m_df[, ditch_presence := fifelse(best == 1, "no",
                                fifelse(best == 2, "yes", NA_character_))]

# Extract results
stat <- as.data.frame(res_m$stat)
lower <- as.data.frame(res_m$lowerCI)
upper <- as.data.frame(res_m$upperCI)

# Compute the rice-wetland difference
assoc <- stat$no  # This is positive toward rice
assoc[is.na(assoc)] <- 0

# Combine into one data frame
df_plot <- data.frame(
  species = rownames(stat),
  assoc = assoc,
  lowerCI = lower$no,
  upperCI = upper$no
)

# Merge with significance info
df_plot <- left_join(df_plot, res_sig_m_df, by = c("species" = "rn"))



# Add direction (rice vs wetland)
df_plot <- df_plot %>%
  mutate(
    assoc_dir = ifelse(assoc > 0, "Absent", "Present"),
    sig = case_when(
      psidak <= 0.001 ~ "***",
      psidak <= 0.01 ~ "**",
      psidak <= 0.05  ~ "*",
      TRUE           ~ "ns"
    ),
    sig_fdr = case_when(
      pfdr <= 0.001 ~ "***",
      pfdr <= 0.01 ~ "**",
      pfdr <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

# Sort by association strength
#df_plot$species <- factor(df_plot$species, levels = df_plot$species[order(df_plot$assoc)])



absent_col <- "gray38"
present_col <- "lightgoldenrod3"

# 1. Reshape into long format
df_long <- df_plot %>%
  dplyr::select(species, assoc, lowerCI, upperCI, assoc_dir, sig, sig_fdr) %>%
  pivot_longer(cols = c(sig, sig_fdr),
               names_to = "correction", values_to = "signif") %>%
  mutate(
    correction = case_when(
      correction == "sig" ~ "Sidak",
      correction == "sig_fdr" ~ "FDR"
    ),
    # x-offset to avoid overlap: -0.015, 0, +0.015
    x_offset = case_when(
      correction == "Sidak" ~ -0.05,
      correction == "FDR"   ~ 0.05
    ),
    assoc_pos = assoc + x_offset
  )

df_long <- df_long %>%
  mutate(sp_group = species_to_group[species])

# Define taxonomic chronological order of groups
group_order <- c("Gastropoda", "Bivalvia", "Amphipod", "Isopoda",
                 "Ephemeroptera", "Odonata", "Megaloptera",
                 "Trichoptera", "Coleoptera", "Diptera", "Hemiptera")


# Reorder species within the plot by group first, then assoc
df_plot <- df_plot %>%
  mutate(group = factor(species_to_group[species], levels = group_order)) %>% # preserve order
  arrange(group, desc(assoc))

# Reset species as factor with the new order
df_plot$species <- factor(df_plot$species, levels = df_plot$species)


# Recompute y positions for vertical group lines
df_plot$y_pos <- as.numeric(factor(df_plot$species))

group_lines <- df_plot %>%
  group_by(group) %>%
  summarize(
    y_min = min(y_pos),
    y_max = max(y_pos),
    y_mid = (min(y_pos) + max(y_pos)) / 2
  )

# Set x position for group labels
label_x <- -1  # adjust left of the plot

##
p_indicator_Present <- ggplot(df_plot, aes(x = assoc, y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.2,
                 color = "gray50", size = 0.8) +
  
  # Add three significance symbols per species
  geom_point(data = df_long,
             aes(x = assoc_pos, y = species,
                 shape = signif, color = correction),
             size = 3.5) +
  
  # Species labels
  geom_text(data = df_plot[df_plot$assoc_dir == "Present", ],
            aes(label = gsub("_", " ", species)),
            hjust = 1.1, vjust = 0.3,
            x = df_plot$lowerCI[df_plot$assoc_dir == "Present"],
            size = 4, fontface = "bold.italic", color = present_col) +
  geom_text(data = df_plot[df_plot$assoc_dir == "Absent", ],
            aes(label = gsub("_", " ", species)),
            hjust = -0.1, vjust = 0.3,
            x = df_plot$upperCI[df_plot$assoc_dir == "Absent"],
            size = 4, fontface = "bold.italic", color = absent_col) +
  
  # Group annotation lines
  geom_segment(data = group_lines,
               aes(x = label_x, xend = label_x, y = y_min, yend = y_max),
               color = "gray30", size = 1) +
  
  # Group labels (no exclamation marks)
  geom_text(data = group_lines,
            aes(x = label_x - 0.05, y = (y_min + y_max)/2, label = group),
            hjust = 1, fontface = "bold", color = "gray30", size = 4)+
  
  # Shapes & colors for significance
  scale_shape_manual(values = c("ns" = 39, "*" = 15, "**" = 17, "***" = 16)) +
  scale_color_manual(values = c("Sidak" = "gold3",
                                "Holm"  = "purple",
                                "FDR"   = "darkgreen")) +
  annotate("segment", x = 0, xend = max(df_plot$assoc) * 0.9, 
           y = length(unique(df_plot$species)) + 4, 
           yend = length(unique(df_plot$species)) + 4, 
           arrow = arrow(length = unit(0.5, "cm")), 
           color = absent_col, size = 3) +
  annotate("text", x = max(df_plot$assoc) * 1, 
           y = length(unique(df_plot$species)) + 4.5, 
           label = "Absent", color = absent_col, hjust = 0, vjust = 0, size = 4.5) +
  
  annotate("segment", x = 0, xend = min(df_plot$assoc) * 0.9, 
           y = length(unique(df_plot$species)) + 4, 
           yend = length(unique(df_plot$species)) + 4, 
           arrow = arrow(length = unit(0.5, "cm")), 
           color = present_col, size = 3) +
  annotate("text", x = min(df_plot$assoc) * 1, 
           y = length(unique(df_plot$species)) + 4.5, 
           label = "Present", color = present_col, hjust = 1, vjust= 0,size = 4.5)+
  
  theme_minimal() +
  labs(
    x = "Point-Biserial Correlation Coefficient",
    y = NULL,
    color = "Correction method",
    shape = "Significance level"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 100, 10, 100, unit = "pt"),
    legend.background = element_rect(fill = "white", colour = "white")
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))

p_indicator_Present


#p_indicator_Present + p_indicator_ditch+plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")&theme(legend.position = "bottom")





#########
# REDUCED (rice CENTER only) with Filtered out side samples
#function to exclude 0

sp_data_r1 <- sp_data_f %>% filter(site_in == "center")
sp_data_r1 <- sp_data_r1 %>%
  dplyr::select(all_of(meta_cols_f), where(~ any(. != 0, na.rm = TRUE)))

sp_data_r <- sp_data_r1[,(30:ncol(sp_data_r1))]
env_r <- sp_data_r1[,(1:29)]
env_r[,c(9,11:13,15)] <- log1p(env_r[,c(9,11:13,15)])
env_r[,c(9,11:15)] <- scale(env_r[,c(9,11:15)])

env_only_r <- env_r[,c(9,11:15)]


#calculate effect size and CI
set.seed(1)
res_r <- strassoc(sp_data_r, sp_data_r1$ditch_present, func="r.g", nboot.ci = 1000)
res_r

# permutation-based p-values (with stratification)
set.seed(1)
res_sig_r <- signassoc(sp_data_r, cluster=sp_data_r1$ditch_present, mode = 1, control = how(
  blocks = sp_data_r1$site_name,  # Only restrict permutations to within-site
  within = Within(type = "free"),
  nperm = 9999))

res_sig_r_df <- as.data.table(res_sig_r, keep.rownames=TRUE)
res_sig_r_df$species <- res_sig_r_df$rn

#corrections
#res_sig_r_df$pholm <- p.adjust(res_sig_r_df$psidak, method = "holm")
res_sig_r_df$pfdr <- p.adjust(res_sig_r_df$psidak, method = "BH")


res_sig_r_df[psidak<=0.05, ]


#Add habitat information automatically
res_sig_r_df[, ditch_presence := fifelse(best == 1, "no",
                                         fifelse(best == 2, "yes", NA_character_))]

# Extract results
stat_r <- as.data.frame(res_r$stat)
lower_r <- as.data.frame(res_r$lowerCI)
upper_r <- as.data.frame(res_r$upperCI)

# Compute the rice-wetland difference
assoc_r <- stat_r$no  # This is positive toward rice
assoc_r[is.na(assoc_r)] <- 0

# Combine into one data frame
df_plot_r <- data.frame(
  species = rownames(stat_r),
  assoc = assoc_r,
  lowerCI = lower_r$no,
  upperCI = upper_r$no
)

# Merge with significance info
df_plot_r <- left_join(df_plot_r, res_sig_r_df, by = c("species" = "rn"))

# Add direction (rice vs wetland)
df_plot_r <- df_plot_r %>%
  mutate(
    assoc_dir = ifelse(assoc > 0, "Absent", "Present"),
    sig = case_when(
      psidak <= 0.001 ~ "***",
      psidak <= 0.01 ~ "**",
      psidak <= 0.05  ~ "*",
      TRUE           ~ "ns"
    ),
    sig_fdr = case_when(
      pfdr <= 0.001 ~ "***",
      pfdr <= 0.01 ~ "**",
      pfdr <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )




absent_col <- "gray38"
present_col <- "gold3"

df_long_r <- df_plot_r %>%
  dplyr::select(species, assoc, lowerCI, upperCI, assoc_dir, sig, sig_fdr) %>%
  pivot_longer(cols = c(sig, sig_fdr),
               names_to = "correction", values_to = "signif") %>%
  mutate(
    correction = case_when(
      correction == "sig" ~ "Sidak",
      correction == "sig_fdr" ~ "FDR"
    ),
    # x-offset to avoid overlap: -0.015, 0, +0.015
    x_offset = case_when(
      correction == "Sidak" ~ -0.05,
      correction == "FDR"   ~ 0.05
    ),
    assoc_pos = assoc + x_offset
  )

df_long_r <- df_long_r %>%
  mutate(sp_group = species_to_group[species])

# Reorder species within the plot by group first, then assoc
df_plot_r <- df_plot_r %>%
  mutate(group = factor(species_to_group[species], levels = group_order)) %>% # preserve order
  arrange(group, desc(assoc))

# Reset species as factor with the new order
df_plot_r$species <- factor(df_plot_r$species, levels = df_plot_r$species)


# Recompute y positions for vertical group lines
df_plot_r$y_pos <- as.numeric(factor(df_plot_r$species))

group_lines_r <- df_plot_r %>%
  group_by(group) %>%
  summarize(
    y_min = min(y_pos),
    y_max = max(y_pos),
    y_mid = (min(y_pos) + max(y_pos)) / 2
  )

# Set x position for group labels
label_x <- -1  # adjust left of the plot



##
p_indicator_Present_C <- ggplot(df_plot_r, aes(x = assoc, y = species)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = lowerCI, xmax = upperCI), height = 0.2,
                 color = "gray50", size = 0.8) +
  
  # Add three significance symbols per species
  geom_point(data = df_long_r,
             aes(x = assoc_pos, y = species,
                 shape = signif, color = correction),
             size = 3.5) +
  
  # Species labels
  geom_text(data = df_plot_r[df_plot_r$assoc_dir == "Present", ],
            aes(label = gsub("_", " ", species)),
            hjust = 1.1, vjust = 0.3,
            x = df_plot_r$lowerCI[df_plot_r$assoc_dir == "Present"],
            size = 4, fontface = "bold.italic", color = present_col) +
  geom_text(data = df_plot_r[df_plot_r$assoc_dir == "Absent", ],
            aes(label = gsub("_", " ", species)),
            hjust = -0.1, vjust = 0.3,
            x = df_plot_r$upperCI[df_plot_r$assoc_dir == "Absent"],
            size = 4, fontface = "bold.italic", color = absent_col) +
  
  # Group annotation lines
  geom_segment(data = group_lines_r,
               aes(x = label_x, xend = label_x, y = y_min, yend = y_max),
               color = "gray30", size = 1) +
  
  # Group labels (no exclamation marks)
  geom_text(data = group_lines_r,
            aes(x = label_x - 0.05, y = (y_min + y_max)/2, label = group),
            hjust = 1, fontface = "bold", color = "gray30", size = 4)+
  
  # Shapes & colors for significance
  scale_shape_manual(values = c("ns" = 39, "*" = 15, "**" = 17, "***" = 16)) +
  scale_color_manual(values = c("Sidak" = "gold3",
                                "Holm"  = "purple",
                                "FDR"   = "darkgreen")) +
  annotate("segment", x = 0, xend = max(df_plot_r$assoc) * 0.9, 
           y = length(unique(df_plot_r$species)) + 4, 
           yend = length(unique(df_plot_r$species)) + 4, 
           arrow = arrow(length = unit(0.5, "cm")), 
           color = absent_col, size = 3) +
  annotate("text", x = max(df_plot_r$assoc) * 1, 
           y = length(unique(df_plot_r$species)) + 4.5, 
           label = "Absent", color = absent_col, hjust = 0, vjust = 0, size = 4.5) +
  
  annotate("segment", x = 0, xend = min(df_plot_r$assoc) * 0.9, 
           y = length(unique(df_plot_r$species)) + 4, 
           yend = length(unique(df_plot_r$species)) + 4, 
           arrow = arrow(length = unit(0.5, "cm")), 
           color = present_col, size = 3) +
  annotate("text", x = min(df_plot_r$assoc) * 1, 
           y = length(unique(df_plot_r$species)) + 4.5, 
           label = "Present", color = present_col, hjust = 1, vjust= 0,size = 4.5)+
  
  theme_minimal() +
  labs(
    x = "Point-Biserial Correlation Coefficient",
    y = NULL,
    color = "Correction method",
    shape = "Significance level"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 100, 10, 100, unit = "pt"),
    legend.background = element_rect(fill = "white", colour = "white")
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))

p_indicator_Present_C



p_indicator_Present + p_indicator_Present_C + p_indicator_ditch+plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")&theme(legend.position = "bottom")




#### Result table with GROUPING:

###
## result table
result_tab_m <- df_plot %>% dplyr::select(assoc_dir, group, species, assoc, psidak, pfdr)%>%arrange(assoc_dir, group, species)
result_tab_r <- df_plot_r %>% dplyr::select(assoc_dir, group, species, assoc, psidak, pfdr)%>%arrange(assoc_dir, group, species)
result_tab_d <- df_plot_d %>% dplyr::select(assoc_dir, group, species, assoc, psidak, pfdr)%>%arrange(assoc_dir, group, species)
