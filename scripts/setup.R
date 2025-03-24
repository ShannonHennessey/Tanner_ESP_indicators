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
## Set years
current_year <- 2024
years <- 1998:current_year


## Pull Tanner specimen data
tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")

## ID corner stations
# Read in station general definitions, grid cell centroid coordinates
stations <- read.csv("Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/lookup_tables/station_lookup.csv")

# Define corner stations
corners <- stations %>% 
           filter(STATION_TYPE == "MTCA_CORNER") %>%
           pull(STATION_ID)




