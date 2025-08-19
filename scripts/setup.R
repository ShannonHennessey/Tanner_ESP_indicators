## Purpose: To serve as a central script to load packages, common variables, 
##          and lookups.


## Load packages
library(crabpack)
library(tidyverse)
library(tidync)
library(lubridate)
library(sf)
library(httr)
library(akgfmaps)
library(rnaturalearth) 


## Set data directory
data_dir <- "Y:/KOD_Research/Hennessey/Tanner_ESP/data/"
fig_dir <- "Y:/KOD_Research/Hennessey/Tanner_ESP/figures/"


## Set years
current_year <- 2025
years <- 1982:current_year


## Pull Tanner specimen data
tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")


## Pull size at 50% probability of terminal molt
mat_size <- get_male_maturity(species = "TANNER", 
                              region = "EBS")$model_parameters %>% 
            dplyr::select(-c("A_EST", "A_SE")) %>%
            rename(MAT_SIZE = B_EST, 
                   STD_ERR = B_SE) %>%
            right_join(., expand_grid(YEAR = 1988:current_year,
                                      SPECIES = "TANNER", 
                                      REGION = "EBS",
                                      DISTRICT = c("ALL", "E166", "W166"))) %>%
            group_by(DISTRICT) %>%
            # assign mean cutline to missing years
            mutate(MAT_SIZE = ifelse(is.na(MAT_SIZE), mean(MAT_SIZE, na.rm = TRUE), MAT_SIZE)) %>%
            arrange(DISTRICT, YEAR)

## ID corner stations
# Read in station general definitions, grid cell centroid coordinates
stations <- read.csv("Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/lookup_tables/station_lookup.csv")

# Define corner stations
corners <- stations %>% 
           filter(STATION_TYPE == "MTCA_CORNER") %>%
           pull(STATION_ID)


## Update Tanner core area
source("./scripts/get_core_area.R")
# - make this into a yearly metric??


