## Purpose: To create a master .csv file of Tanner crab ecosystem indicators for 
##          R Markdown


# ## Load packages
library(tidyverse)
# # library(corrplot)
# # library(cowplot)
# # library(mgcv)


## Set data directory
data_dir <- "Y:/KOD_Research/Hennessey/Tanner_ESP/data/"

## Read in indicator time series
# Chlorophyll-a concentration
chla <- read.csv("./outputs/chla_concentration.csv") # Contributor: Matt Callahan

# Wind -- northeasterly favors larval advection...
# use cross-shelf component only? Negative would be "northeasterly"
wind <- read.csv("./outputs/wind_NCEP.csv") # Contributor: Tyler Hennon

# Mean summer bottom temperature, summer cold pool extent, 
# mean summer surface temperature, Aleutian Low
env <- read.csv("./outputs/temp_coldpool_AL.csv")

# North Pacific Index
npi <- read.csv(paste0(data_dir, "NPI.csv"))

# Sea ice extent
ice <- read.csv(paste0(data_dir, "ice_extent.csv"))

# Summer juvenile Tanner temperature occupancy
occ <- read.csv("./outputs/tanner_temp_occupied.csv")

# # DFA - temperature
# dfa_temp <- read.csv("./outputs/juv_temp_dfa_trend.csv") %>%
#             select(year, dfa_temp)

# Predator density
pred <- read.csv("./outputs/benthic_predator_density.csv")

# Daily summer pacific cod consumption of Tanner crab 
pcod <- read.csv(paste0(data_dir, "pcod_consumption.csv")) # Contributor: Kerim Aydin

# Summer juvenile Tanner disease prevalence
bcd <- read.csv("./outputs/bcd_prevalence.csv") %>%
       dplyr::select(YEAR, SM_PREVALENCE) %>%
       rename(year = YEAR,
              bcd_prevalence = SM_PREVALENCE)

# -- Juvenile cohort progression? or recruitment propagation index?


# Summer benthic invertebrate density
prey <- read.csv("./outputs/benthic_invert_density.csv") %>%
        mutate(total_invert = ifelse(year >= 1988, total_invert, NA))

# Female size at maturity
female_sam <- read.csv("./outputs/female_SAM.csv")
  
# Female clutch fullness
female_clutch <- read.csv("./outputs/clutch_empty.csv")

# Annual Tanner size at terminal molt
male_sam <- read.csv("./outputs/male_term_molt.csv")  # Contributor: Jon Richar
  
# Summer mature male Tanner area occupied (D95)
d95 <- read.csv("./outputs/tanner_area_occupied.csv") %>%
       dplyr::select(YEAR, mature_male) %>%
       rename(year = YEAR, 
              matmale_d95 = mature_male)

# Summer mature male tanner centroid of abundance
cod <- read.csv("./outputs/tanner_centroid_abundance.csv") %>%
       dplyr::select(YEAR, LAT_COD_mature_male, LON_COD_mature_male) %>%
       rename(year = YEAR,
              matmale_cod_lat = LAT_COD_mature_male,
              matmale_cod_lon = LON_COD_mature_male)

# Fraction of stock in BB
frac_bb <- read.csv("./outputs/fraction_bb.csv")

# -- DFA - male spatial? (^^above 3)
# -- Snow/Tanner spatial overlap?



## Combine indices and save output
indicators <- list(chla, env, npi, wind, ice, occ, #dfa_temp, 
                   pred, pcod, bcd, prey, female_sam,
                   female_clutch, male_sam, d95, cod, frac_bb) %>%
              purrr::reduce(., merge, by = c('year'), all = T) %>%
              filter(year >= 1975) %>%
              write_csv("./data/BAS_indicators.csv")



