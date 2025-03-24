## Purpose: To create a master .csv file of Tanner crab ecosystem(?) indicators for 
##          R Markdown


# ## Load packages
# library(tidyverse)
# # library(corrplot)
# # library(cowplot)
# # library(mgcv)


## Read in indicator time series
# Chlorophyll-a concentration
# chla <- read_csv("./outputs/chla.csv") # Contributor: Matt Callahan

# Wind?
# Aleutian Low?
# Summer surface temperature?


# Mean summer bottom temperature, summer cold pool extent
env <- read_csv("./outputs/temp_coldpool.csv")

# Summer juvenile Tanner temperature occupancy
occ <- read_csv("./outputs/tanner_temp_occupied.csv")

# DFA - temperature? (^^above 3)

# # Predator density?
# pred <- read_csv("./outputs/benthic_predator_density.csv")

# Daily summer pacific cod consumption of Tanner crab 
pcod <- read_csv("./data/EBS_codcrab_sum.csv", show_col_types = FALSE) # Contributor: Kerim Aydin

# DFA - predation?


# Summer juvenile Tanner disease prevalence
bcs <- read_csv("./outputs/bcs_prevalence.csv")

# Juvenile cohort progression? or recruitment propagation index?


# # Summer benthic invertebrate density
# invert <- read_csv("./outputs/benthic_invert_density.csv")

# Female size at maturity?
# Female clutch fullness?

# Annual Tanner size at terminal molt
size <- read_csv("./outputs/male_term_molt.csv") # Contributor: Jon Richar

# Summer mature male Tanner area occupied (D95)
d95 <- read_csv("./outputs/tanner_area_occupied.csv")

# Summer mature male tanner centroid of abundance
cod <- read_csv("./outputs/tanner_centroid_abundance.csv")

# Fraction of stock in BB?

# DFA - male spatial? (^^above 3)


# Snow/Tanner spatial overlap?







# combine indices and save output
eco_ind <- invert %>%
           select(YEAR, Total_Benthic) %>%
           rename(Summer_Benthic_Invertebrate_Density_SEBS_Tanner_Survey = Total_Benthic) %>%
           full_join(env %>%
                       select(YEAR, cp_extent, summer_bt) %>%
                       rename(Summer_Cold_Pool_SEBS_Tanner_Survey=cp_extent, Summer_Temperature_Bottom_Tanner_Survey=summer_bt)) %>%
           full_join(d95 %>%
                       select(YEAR, mature_male) %>%
                       rename(Summer_Tanner_Male_Area_Occupied_SEBS_Survey=mature_male)) %>%
           full_join(bcs %>%
                       select(YEAR, Immature) %>%
                       rename(Summer_Tanner_Juvenile_Disease_Prevalence=Immature)) %>%
           full_join(cod %>%
                       select(YEAR, Lon_COD_mature_male) %>%
                       rename(Summer_Tanner_Male_Center_Distribution_SEBS_Survey = Lon_COD_mature_male)) %>%
           full_join(occ %>%
                       select(YEAR, Immature) %>%
                       rename(Summer_Tanner_Juvenile_Temperature_Occupancy = Immature)) %>%
           full_join(size %>%
                       select(Year, SAM_pop) %>%
                       rename(Annual_Tanner_Size_Terminal_Molt_Model = SAM_pop,
                           YEAR = Year)) %>%
           full_join(chla) %>%
           full_join(pcod %>%
                       select(YEAR, Pcod) %>%
                       rename(Summer_Pacific_Cod_Density_Tanner_Survey = Pcod)) %>%
           full_join(pcod %>%
                       select(YEAR, bairdi_extrap) %>%
                       rename(Summer_Pacific_Cod_Consumption_Juvenile_Tanner = bairdi_extrap)) %>%
           rename(year = YEAR) %>%
           filter(year >= 1982) %>%
           arrange(year)
write_csv(eco_ind, "./data/tanner_eco_indicators.csv")



