# Purpose: Pull zero-filled CPUE tables from the RACE Oracle Schema
#          for calculating mean CPUE of groundfish and benthic 
#          invertebrate guilds 

library(gapindex)

# Connect to AFSC Oracle Database
channel <- gapindex::get_connected()

# Manipulate with SQL using AKFIN data before bringing into R ------------------
  # This will pull a final, completely formatted table for only EBS and NBS
  # FYI: This is a long run time (~4.5 hours), and you must be connected to VPN!
start_time <- Sys.time()
a <- RODBC::sqlQuery(channel = channel, # NOT RACEBASE.HAUL
                     query = paste0("SELECT
                                    cr.CRUISEJOIN,
                                    cr.CRUISE,
                                    cr.YEAR,
                                    cr.SURVEY_DEFINITION_ID,
                                    cr.SURVEY_NAME,
                                    cr.VESSEL_ID,
                                    cr.VESSEL_NAME,
                                    cp.HAULJOIN,
                                    cp.SPECIES_CODE,
                                    tt.SPECIES_NAME,
                                    tt.COMMON_NAME,
                                    cp.WEIGHT_KG,
                                    cp.COUNT,
                                    cp.AREA_SWEPT_KM2,
                                    cp.CPUE_KGKM2,
                                    cp.CPUE_NOKM2,
                                    -- cp.CPUE_KGKM2/100 AS WTCPUE,
                                    -- cp.CPUE_NOKM2/100 AS NUMCPUE,
                                    hh.HAUL,
                                    hh.STATION,
                                    hh.LATITUDE_DD_START,
                                    hh.LATITUDE_DD_END,
                                    hh.LONGITUDE_DD_START,
                                    hh.LONGITUDE_DD_END
                                    FROM GAP_PRODUCTS.AKFIN_HAUL hh
                                    LEFT JOIN GAP_PRODUCTS.AKFIN_CRUISE cr
                                    ON hh.CRUISEJOIN = cr.CRUISEJOIN
                                    LEFT JOIN GAP_PRODUCTS.AKFIN_CPUE cp
                                    ON hh.HAULJOIN = cp.HAULJOIN
                                    LEFT JOIN GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION tt
                                    ON cp.SPECIES_CODE = tt.SPECIES_CODE
                                    WHERE SURVEY_DEFINITION_ID IN (143, 98) -- 143 NBS, 98 EBS
                                    AND tt.SURVEY_SPECIES = 1;")) 
end_time <- Sys.time()
cat(paste("Time Elapsed:", round(end_time - start_time, 2), 
          units(end_time - start_time), "\n\n"))

write.csv(a, "./data/gf_cpue_timeseries.csv")

# Alternatively, you an just download the files and manipulate locally ---------
# this will pull all standard RACE survey data (e.g., also GOA, AI, BSS)

### Load data files from Oracle 
locations <- c("GAP_PRODUCTS.AKFIN_HAUL", 
               "GAP_PRODUCTS.AKFIN_CRUISE", 
               "GAP_PRODUCTS.AKFIN_CPUE", 
               "GAP_PRODUCTS.TAXONOMIC_CLASSIFICATION")

for(i in 1:length(locations)){
  print(locations[i])
  a <- RODBC::sqlQuery(channel = channel, # NOT RACEBASE.HAUL
                       query = paste0("SELECT * FROM ", locations[i], "; ")) 
  write.csv(a, paste0("./data/", tolower(gsub(pattern = '.', 
                                              replacement = "_", 
                                              x = locations[i], 
                                              fixed = TRUE)), ".csv"),
            row.names = FALSE)
}

### Load data from local data folder 
print("Load oracle data")
a <- list.files(path = here::here("data"), 
                full.names = TRUE, recursive = FALSE, pattern = "gap_products")

for (i in 1:length(a)){
  b <- read_csv(file = a[i], show_col_types = FALSE)
  temp <- strsplit(x = a[i], split = "/")
  temp <- toupper(gsub(pattern = "\\.csv", replacement = "", x = temp[[1]][length(temp[[1]])]))
  assign(x = temp, value = b)
}

### Wrangle code in R 
dat <- dplyr::left_join(GAP_PRODUCTS_AKFIN_HAUL, 
                        GAP_PRODUCTS_AKFIN_CRUISE) %>% 
       dplyr::left_join(GAP_PRODUCTS_AKFIN_CPUE) %>% 
       dplyr::left_join(GAP_PRODUCTS_TAXONOMIC_CLASSIFICATION %>% filter(SURVEY_SPECIES == 1)) %>% 
       dplyr::filter(SURVEY_DEFINITION_ID %in% c(143, 98)) %>% 
       dplyr::select(CRUISEJOIN, CRUISE, YEAR, SURVEY_DEFINITION_ID, 
                     SURVEY_NAME, VESSEL_ID, VESSEL_NAME, HAULJOIN, 
                     SPECIES_CODE, SPECIES_NAME, COMMON_NAME, WEIGHT_KG, 
                     COUNT, AREA_SWEPT_KM2, CPUE_KGKM2, CPUE_NOKM2,
                     HAUL, STATION) %>% 
       dplyr::mutate(WTCPUE = CPUE_KGKM2/100,
                     NUMCPUE = CPUE_NOKM2/100)

write.csv(dat, "./data/gf_cpue_timeseries.csv", row.names = FALSE)

test <- read.csv("./data/gf_cpue_timeseries_2024.csv")

# Manipulate with SQL using FOSS data before bringing into R -------------------
# This will pull a final, completely formatted table for only EBS and NBS
start_time <- Sys.time()
dat2 <- RODBC::sqlQuery(channel = channel, # NOT RACEBASE.HAUL
                       query = paste0("SELECT
                                       hh.CRUISEJOIN,
                                       hh.CRUISE,
                                       hh.YEAR,
                                       --hh.SURVEY_DEFINITION_ID,
                                       hh.SURVEY,
                                       hh.VESSEL_ID,
                                       hh.VESSEL_NAME,
                                       hh.HAULJOIN,
                                       hh.HAUL,
                                       hh.STATION,
                                       hh.LATITUDE_DD_START,
                                       hh.LATITUDE_DD_END,
                                       hh.LONGITUDE_DD_START,
                                       hh.LONGITUDE_DD_END,
                                       hh.AREA_SWEPT_KM2,
                                       cc.SPECIES_CODE,
                                       tt.SCIENTIFIC_NAME,
                                       tt.COMMON_NAME,
                                       cc.WEIGHT_KG,
                                       cc.COUNT,
                                       cc.CPUE_KGKM2,
                                       cc.CPUE_NOKM2
                                       FROM GAP_PRODUCTS.FOSS_HAUL hh
                                       LEFT JOIN GAP_PRODUCTS.FOSS_CATCH cc
                                       ON hh.HAULJOIN = cc.HAULJOIN
                                       LEFT JOIN GAP_PRODUCTS.FOSS_SPECIES tt
                                       ON cc.SPECIES_CODE = tt.SPECIES_CODE
                                       WHERE SURVEY_DEFINITION_ID IN (143, 98);")) # 143 NBS, 98 EBS
end_time <- Sys.time()
cat(paste("Time Elapsed:", round(end_time - start_time, 2), 
          units(end_time - start_time), "\n\n"))

write.csv(dat, "./data/gf_cpue_timeseries.csv", row.names = FALSE)


