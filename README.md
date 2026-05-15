    STUDY DESCRIPTION
    This dataset originates from the study "Rice paddies promote diverse and distinct aquatic invertebrate communities in agroecosystems", which investigates the role of pesticide-free rice paddies as aquatic habitats in the Swiss lowlands.

    The study compares aquatic macroinvertebrate communities and environmental conditions between rice paddies and nearby natural wetlands over two years (2022–2023), with two sampling periods per year (July and August).

    A total of 11 sites were analysed, each consisting of a rice paddy and a nearby wetland.

    Rice paddies represent temporary aquatic systems that are flooded during the growing season (spring–summer) and drained after harvest, whereas wetlands vary in hydroperiod, structure, and vegetation.


    STUDY DESIGN
    Sampling design:
    - two years: 2022–2023
    - two sampling rounds per year: July, August
    - paired design: rice paddy vs. wetland per site
    - additional structure: presence/absence of side ditches in rice paddies

    Macroinvertebrates were sampled using standardized net sampling and counts were converted to densities (individuals per m²) to allow comparison across sites.

    Environmental variables were measured in parallel (e.g. conductivity, oxygen, pH, temperature, depth).

    RELATION TO R ANALYSIS SCRIPTS
    The dataset is designed to be used together with the provided R scripts, which follow the manuscript structure:

    1_Environment.R
    1_2_Environment_ricefields.R
        -> analysis of abiotic/environmental variables

    2_diversity.R
    2_2_diversity_rice_fields_.R
        -> diversity analyses (alpha, gamma, density, Shannon)

    3_2_ISA_rice_fields_.R
    3_ISA_dbRDA_wetland_rice_plot.R
        -> community composition analysis
           - indicator species analysis (ISA)
           - multivariate ordination (dbRDA)

    Scripts are numbered according to manuscript sections.

        DATASET STRUCTURE

    Columns are organised in two main blocks:

    1) site_name → water_stop
       = metadata + environmental + management variables

    2) Coenangrionidae → last column
       = taxa abundance data (counts per sample)


    COLUMN GROUPS in Excel file: 2022_2023_macroinvRADICAL_abio_env_20250407.xlxs

    [METADATA]
    site_name              : site identifier
    env                    : environment type (e.g. rice, wetland)
    pond_#                 : pond identifier within site
    sampling_round         : sampling event index
    year                   : sampling year (-)
    date                   : sampling date (YYYY-MM-DD or DD/MM/YYYY)
    site_in                : sampling position (e.g. centre, side)
    total_subsamples       : number of subsamples (count)
    ID by                  : person who identified sample


    [COMMUNITY SUMMARY]
    Sobs                   : observed taxa richness (count)
    n                      : total individuals (count)


    [GEOGRAPHY]
    canton                 : Swiss canton
    ch_NS                  : region (North / South)
    area_ha                : area (ha)
    distance_rice_wetland_m: distance between habitats (m)


    [ABIOTIC VARIABLES]
    cond                   : conductivity (µS/cm)
    depth                  : water depth (m)
    NO                     : nitrogen concentration (mg/L)
    oxy                    : dissolved oxygen (mg/L)
    ph                     : pH (-)
    PO                     : phosphorus concentration (mg/L)
    temp                   : water temperature (°C)


    [AGRICULTURAL MANAGEMENT]
    N_fertiliser_name      : fertiliser name
    fertiliser_type        : fertiliser type (e.g. mineral, organic)
    total_N_per_ha         : nitrogen input (kg N/ha)
    fertiliser_latest_date : last fertilisation date

    soil_prep_type         : soil preparation method
    soil_prep_date         : soil preparation date
    soil_prep_depth_cm     : soil preparation depth (cm)

    rice_transplant        : transplanting date
    rice_harvest           : harvest date


    [WATER MANAGEMENT]
    water_overwinter       : overwintering water presence (yes/no)
    water_start            : flooding start date
    water_stop             : flooding end date


    [TAXA DATA]
    Columns from "Coenangrionidae" to the final column.
    - each column = one taxon (family / genus / species)
    - each value  = count of individuals (ind.)
    - dataset is in wide format (samples = rows, taxa = columns)

