log_cpue <- function(crab_data = NULL,
                      species = NULL,
                      region = c("EBS", "NBS")[1],
                      district = c("ALL", "BB", "NORTH", "NS", "PRIB", "STMATT", "UNSTRAT", "E166", "W166", "166TO173")[1],
                      years = NULL,
                      sex = NULL,
                      size_min = NULL,
                      size_max = NULL,
                      crab_category = NULL,
                      female_maturity = c("morphometric", "cutline")[1],
                      shell_condition = NULL,
                      egg_condition = NULL,
                      clutch_size = NULL,
                      bin_1mm = FALSE,
                      replace_retow = TRUE,
                      # rm_corners = FALSE,
                      output = c("cpue", "bioabund")[1]){
  
  
  ## Check inputs
  if(is.null(x = crab_data)){
    stop(paste0("Must provide argument `crab_data` a named list from",
                " `crabpack::get_specimen_data()`."))
  }
  
  
  # Can only handle 1 species at a time
  if(length(species) > 1){
    stop("Argument `species` must be of length = 1.")
  }
  
  
  ## ERROR: check species in data == species called in function
  if(unique(crab_data$specimen$SPECIES) != species){
    stop("Argument `species` does not match species in crab_data.")
  }
  
  
  ## ERROR: check that district(s) specified are correct for given species
  if(TRUE %in% (species %in% c("RKC", "HAIR")) &
     TRUE %in% (!district %in% c("ALL", "BB", "PRIB", "NORTH")) &
     region == "EBS"){
    stop(paste0("Valid `district` options for EBS ", species, " are:",
                " c('ALL', 'BB', 'PRIB', 'NORTH').",
                " Please set `district` to one or many of these options."))
  }
  
  if(species == "BKC" &
     TRUE %in% (!district %in% c("ALL", "PRIB", "STMATT", "UNSTRAT")) &
     region == "EBS"){
    stop(paste0("Valid `district` options for EBS ", species, " are:",
                " c('ALL', 'PRIB', 'STMATT', 'UNSTRAT').",
                " Please set `district` to one or many of these options."))
  }
  
  if(species == "TANNER" & TRUE %in% (!district %in% c("ALL", "E166", "W166", "166TO173")) &
     region == "EBS"){
    stop(paste0("Valid `district` options for EBS ", species, " are:",
                " c('ALL', 'E166', 'W166', '166TO173').",
                " Please set `district` to one or many of these options."))
  }
  
  if(species %in% c("SNOW", "HYBRID") & TRUE %in% (!district %in% c("ALL")) &
     region == "EBS"){
    stop(paste0("Valid `district` options for EBS ", species, " are:",
                " c('ALL'). Please set `district` to this option."))
  }
  
  
  if(TRUE %in% (species %in% c("RKC")) &
     TRUE %in% (!district %in% c("ALL", "NS")) &
     region == "NBS"){
    stop(paste0("Valid `district` options for NBS ", species, " are:",
                " c('ALL', 'NS').",
                " Please set `district` to one or many of these options."))
  }
  
  if(TRUE %in% (species %in% c("BKC", "TANNER", "SNOW", "HYBRID", "HAIR")) &
     TRUE %in% (!district %in% c("ALL")) &
     region == "NBS"){
    stop(paste0("Valid `district` options for NBS ", species, " is:",
                " c('ALL').",
                " Please set `district` to this option."))
  }
  
  ## If "ALL" is specified in `district`, just set district to all
  if("ALL" %in% district){
    district <- "ALL"
  }
  
  
  ## ERROR: Can't calculate both W166 AND 166TO173 for Tanner, they overlap and things may be wonky. Do separately
  if("166TO173" %in% district & length(district > 1)){
    stop(paste0("The 166TO173W district is a special case area for EBS Tanner Crab, and",
                " should not be calculated in conjunction with any other districts due to",
                " spatial overlap with W166. Please set `district` to just `166TO173`",
                " if you wish to query that district. Data may also have to be pulled using",
                " `crabpack::get_specimen_data` specifically for this species and district."))
  }
  
  
  ## ERROR: must specify just one output if length(output > 1)
  if(length(output) > 1){
    stop("Argument `output` must be of length = 1.")
  }
  
  
  ## WARNING: specifying 1mm bins across multiple categories may exceed memory limits
  if(!is.null(crab_category)){
    if(bin_1mm == TRUE & TRUE %in% c(length(crab_category) > 2 | "all_categories" %in% crab_category)){
      
      warning_yn = function(w) {
        cat("Warning message:", '\n', w, '\n')
        cont <- readline('Do you wish to continue? [Y/N] ')
        if(!cont %in% c('Y', 'y')) stop('Aborted by user', call. = FALSE)
      }
      
      warning_yn(paste0("Calculating across multiple categories by 1mm bin may exceed",
                        " R memory limits, especially if metrics are being calculated",
                        " across the entire region. We suggest calculating by individual",
                        " categories instead, and appending the outputs."))
    }
  }
  
  
  
  ## Set variables, define and filter specimen modifiers
  vars <- crabpack::set_variables(crab_data = crab_data,
                                  species = species,
                                  sex = sex,
                                  size_min = size_min,
                                  size_max = size_max,
                                  crab_category = crab_category,
                                  female_maturity = female_maturity,
                                  shell_condition = shell_condition,
                                  egg_condition = egg_condition,
                                  clutch_size = clutch_size,
                                  bin_1mm = bin_1mm)
  
  
  # Pull column names and expand_grid() categories to track throughout calculations
  group_cols <- vars$group_cols
  list2env(vars$expand_combos, .GlobalEnv)
  
  
  # define year if not specified
  current_year <- as.numeric(format(Sys.Date(), "%Y"))
  
  if(is.null(years)){
    years <- c(1975:current_year)
  }
  
  
  # Pull specimen dataframe
  specimen_dat <- vars$specimen_data
  
  
  # Filter by species, make sure species is the same as the data pulled
  specimen_dat <- specimen_dat %>%
    dplyr::filter(SPECIES %in% species)
  if(nrow(specimen_dat) == 0){
    stop(paste0("Selected `species` is not present in the specimen data. Verify",
                " that the species specified in `crabpack::get_specimen_data()`",
                " is the same species specified in `crabpack::calc_cpue()`."))
  }
  
  
  # Specify stock stations
  stock_stations <- crab_data$haul %>%
    # Filter years, region, district
    dplyr::filter(YEAR %in% years,
                  REGION %in% region,
                  !is.na(TOTAL_AREA)) %>%
    {if(!"ALL" %in% district) dplyr::filter(., DISTRICT %in% district) else .} %>%
    # make dummy HT for retow station tracking, all HT = 3 except 17
    dplyr::mutate(HT = ifelse(HAUL_TYPE == 17, 17, 3),
                  LATITUDE = MID_LATITUDE,
                  LONGITUDE = MID_LONGITUDE) %>%
    dplyr::select(HAULJOIN, REGION, YEAR, STATION_ID, HT, AREA_SWEPT,
                  LATITUDE, LONGITUDE, DISTRICT, STRATUM, TOTAL_AREA) %>%
    dplyr::distinct()
  
  # # If want to remove MTCA corner stations, remove stations and change MTCA stratum to Single
  # if(rm_corners == TRUE){
  #   # remove corners within MTCA strata
  #   stock_stations <- stock_stations %>%
  #                     dplyr::filter(!(STRATUM %in% c("PRIB_MTCA", "STMATT_MTCA", "STMATT_MTCA_E") & nchar(STATION_ID) == 6))
  #
  #   # integrate MTCA stratum into the encompassing stratum, based on species
  #   if(species == "BKC"){
  #     stock_stations <- stock_stations %>%
  #                       dplyr::mutate(STRATUM = case_when(STRATUM == "PRIB_MTCA" ~ "PRIB",
  #                                                         STRATUM == "STMATT_MTCA" ~ "STMATT",
  #                                                         TRUE ~ STRATUM))
  #   }
  #
  #   if(species %in% c("RKC", "HAIR")){
  #     # STMATT_MTCA corners --> NORTH for these species, not sure how to remove them
  #     # without removing other accidental corners in older part of ts. Could manually specify
  #     # but not sure if anyone wants/needs these data.
  #     stock_stations <- stock_stations %>%
  #                       dplyr::mutate(STRATUM = dplyr::case_when(STRATUM == "PRIB_MTCA" ~ "PRIB",
  #                                                                TRUE ~ STRATUM))
  #
  #     warning(paste0("Corner stations within the St. Matthew district are considered",
  #                    " part of the Northern district for this species and have not",
  #                    " been removed. Please submit an issue at",
  #                    " github.com/AFSC-Shellfish-Assessment-Program/crabpack/issues",
  #                    " if this data modification is needed."))
  #
  #   }
  #
  #   if(species == "TANNER"){
  #     if(district == "166TO173"){
  #       stock_stations <- stock_stations %>%
  #                         dplyr::mutate(STRATUM = dplyr::case_when(STRATUM == "PRIB_MTCA" ~ "166TO173",
  #                                                                  STRATUM == "STMATT_MTCA_E" ~ "166TO173",
  #                                                                  TRUE ~ STRATUM))
  #     } else{
  #       stock_stations <- stock_stations %>%
  #                         dplyr::mutate(STRATUM = dplyr::case_when(STRATUM == "PRIB_MTCA" ~ "W166",
  #                                                                  STRATUM == "STMATT_MTCA" ~ "W166",
  #                                                                  TRUE ~ STRATUM))
  #     }
  #   }
  #
  #   if(species %in% c("SNOW", "HYBRID")){
  #     stock_stations <- stock_stations %>%
  #                       dplyr::mutate(STRATUM = dplyr::case_when(STRATUM == "PRIB_MTCA" ~ "EBS_SINGLE",
  #                                                                STRATUM == "STMATT_MTCA" ~ "EBS_SINGLE",
  #                                                                TRUE ~ STRATUM))
  #   }
  # }
  
  ## Calculate CPUE by STATION, YEAR, and other biometrics ---------------------
  cpue <- specimen_dat %>%
    # make dummy HT for retow station tracking, all HT = 3 except 17
    dplyr::mutate(HT = ifelse(HAUL_TYPE == 17, 17, 3)) %>%
    # filter hauls in relevant stock stations
    dplyr::filter(HAULJOIN %in% stock_stations$HAULJOIN) %>%
    dplyr::mutate(COUNT = SAMPLING_FACTOR,
                  CPUE = SAMPLING_FACTOR/AREA_SWEPT,
                  CPUE_MT = (SAMPLING_FACTOR * CALCULATED_WEIGHT) / AREA_SWEPT / 1000 / 1000,
                  CPUE_LBS = (SAMPLING_FACTOR * CALCULATED_WEIGHT) * 0.0022046226218488 / AREA_SWEPT) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'HT', 'REGION', 'STATION_ID', 'SEX_TEXT', group_cols)))) %>%
    dplyr::summarise(COUNT = sum(COUNT),
                     CPUE = sum(CPUE),
                     CPUE_MT = sum(CPUE_MT),
                     CPUE_LBS = sum(CPUE_LBS))
  
  
  ## Join to zero catch stations, summarize ------------------------------------
  station_cpue <- cpue %>%
    # dplyr::group_by(YEAR) %>%
    dplyr::right_join(tidyr::expand_grid(SEX_TEXT = sex_combos,
                                         CATEGORY = category_combos,
                                         SHELL_TEXT = shell_combos,
                                         EGG_CONDITION_TEXT = egg_combos,
                                         CLUTCH_TEXT = clutch_combos,
                                         SIZE_1MM = bin_combos,
                                         HT = unique(stock_stations$HT),
                                         stock_stations %>%
                                           dplyr::select(YEAR, REGION, STATION_ID,
                                                         DISTRICT, STRATUM, TOTAL_AREA)),
                      by = c('YEAR', 'HT', 'REGION', 'STATION_ID', 'SEX_TEXT', group_cols)) %>%
    tidyr::replace_na(list(COUNT = 0, CPUE = 0, CPUE_MT = 0, CPUE_LBS = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(stock_stations, ., by = c('REGION', 'YEAR', 'STATION_ID', 'HT', 'DISTRICT', 'STRATUM', 'TOTAL_AREA')) %>%
    dplyr::select(dplyr::all_of(c("YEAR", "HT", "STATION_ID", "SEX_TEXT", group_cols,
                                  "COUNT", "CPUE", "CPUE_MT", "CPUE_LBS",
                                  "REGION", "DISTRICT", "STRATUM", "TOTAL_AREA"))) %>%
    distinct()
  
  
  # add sex_text to group_cols if still needed (ie. specified in function)
  if(!is.null(sex)){
    group_cols <- append(group_cols, "SEX_TEXT")
  }
  
  # remove `sex_combos` if `sex = NULL`
  if(is.null(sex)){
    sex_combos <- NA
  }
  
  # REMOVE NONSENSICAL CRAB CATEGORY COMBOS REMNANT FROM EXPAND_GRID......
  # if we end up making egg_condition and clutch_size not dependent on female only,
  # will probably have to remove those too
  if(!is.null(crab_category)){
    station_cpue <- station_cpue %>%
      dplyr::filter(!((CATEGORY %in% c("mature_male", "large_male", "immature_male", "small_male",
                                       "legal_male", "sublegal_male", "preferred_male") & SEX_TEXT == "female") |
                        (CATEGORY %in% c("mature_female", "immature_female", "female") & SEX_TEXT == "male")))
  }
  
  
  # Filter out STATION E-11 in year 2000 for BBRKC males because it wasn't sampled in leg 1
  # Keep entry in with NA values
  if(species == "RKC" & region == "EBS" & TRUE %in% (district %in% c("ALL", "BB"))){
    station_cpue <- station_cpue %>%
      dplyr::mutate(COUNT = ifelse((YEAR == 2000 & STATION_ID == "E-11" & SEX_TEXT == "male"), NA, COUNT),
                    CPUE = ifelse((YEAR == 2000 & STATION_ID == "E-11" & SEX_TEXT == "male"), NA, CPUE),
                    CPUE_MT = ifelse((YEAR == 2000 & STATION_ID == "E-11" & SEX_TEXT == "male"), NA, CPUE_MT),
                    CPUE_LBS = ifelse((YEAR == 2000 & STATION_ID == "E-11" & SEX_TEXT == "male"), NA, CPUE_LBS))
  } else{
    station_cpue <- station_cpue
  }
  
  
  # Replace retow BBRKC --------------------------------------------------------
  if(species == "RKC" & replace_retow != FALSE){
    # Specify retow stations and years for BBRKC
    ## Retow years: 1999 2000 2006 2007 2008 2009 2010 2011 2012 2017 2021
    retows <- stock_stations %>%
      dplyr::filter(HT == 17) %>%
      dplyr::mutate(RETOW = "yes") %>%
      dplyr::select(-c('HAULJOIN', 'AREA_SWEPT', 'LATITUDE', 'LONGITUDE', 'HT'))
    
    # replace female BBRKC with female data from station with HT 17
    station_cpue <- station_cpue %>%
      dplyr::left_join(., retows, by = c('YEAR', 'STATION_ID', 'REGION', 'DISTRICT', 'STRATUM', 'TOTAL_AREA')) %>%
      dplyr::mutate(REMOVE = dplyr::case_when((SEX_TEXT == 'female' & RETOW == "yes" & HT == 3) ~ "remove",
                                              (SEX_TEXT == 'male' & RETOW == "yes" & HT == 17) ~ "remove",
                                              TRUE ~ "keep")) %>%
      dplyr::filter(REMOVE == "keep") %>%
      {if(is.null(sex)) dplyr::mutate(., HT = ifelse(HT == 17, 3, HT)) else .} %>% # change HT to 3 if 17 so male and female data group properly
      dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'STATION_ID', 'HT', group_cols, # can't have unique HT in if want total crab
                                                    'REGION', 'DISTRICT', 'STRATUM', 'TOTAL_AREA')))) %>%
      dplyr::summarise(COUNT = sum(COUNT),
                       CPUE = sum(CPUE),
                       CPUE_MT = sum(CPUE_MT),
                       CPUE_LBS = sum(CPUE_LBS)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(HT = ifelse((YEAR == 2000 & STATION_ID == "E-11"), 17, HT)) #%>% # change HT to 17 so there's a lat/lon joined down the line
    # {if("male" %in% sex_combos) dplyr::add_row(., YEAR = 2000, STATION_ID = "E-11", HT = 17, SEX_TEXT = "male", REGION = "EBS",
    #                                               DISTRICT = "BB", STRATUM = "BB", TOTAL_AREA = 54536,
    #                                               COUNT = NA, CPUE = NA, CPUE_MT = NA, CPUE_LBS = NA)} # if by sex, add NA row for special case missed tow
  } else{
    # remove all HT 17 data and re-summarise to ignore SEX
    station_cpue <- station_cpue %>%
      dplyr::filter(HT != 17) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'STATION_ID', 'HT', group_cols, # keeping HT in for output verification
                                                    'REGION', 'DISTRICT', 'STRATUM', 'TOTAL_AREA')))) %>%
      dplyr::summarise(COUNT = sum(COUNT),
                       CPUE = sum(CPUE),
                       CPUE_MT = sum(CPUE_MT),
                       CPUE_LBS = sum(CPUE_LBS)) %>%
      dplyr::ungroup() #%>%
    # dplyr::add_row(YEAR = 2000, STATION_ID = "E-11", HT = 17, REGION = "EBS", # add NA row for special case missed tow
    #                DISTRICT = "BB", STRATUM = "BB", TOTAL_AREA = 54536,
    #                COUNT = NA, CPUE = NA, CPUE_MT = NA, CPUE_LBS = NA) ## THESE MIGHT HAVE TO BE 0 instead!!!
  }
  
  
  ## CPUE output ---------------------------------------------------------------
  # make final list of groupings to carry over to `calc_bioabund()`
  if(is.null(sex)){
    groups_out <- group_cols[!group_cols %in% c("SEX_TEXT")]
  } else{
    groups_out <- group_cols
  }
  
  # Format output df
  cpue_out <- station_cpue %>%
    dplyr::left_join(., stock_stations, by = c('YEAR', 'STATION_ID', 'HT', # note that even with BB retow, this will give lat/lon for the original haul
                                               'REGION', 'DISTRICT', 'STRATUM', 'TOTAL_AREA')) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SPECIES = species,
                  CPUE = log(CPUE)) %>%
    dplyr::select(dplyr::all_of(c('SPECIES', 'YEAR', 'REGION', 'STATION_ID', 'LATITUDE', 'LONGITUDE',
                                  groups_out, 'DISTRICT', 'STRATUM', 'TOTAL_AREA',
                                  "COUNT", "CPUE", "CPUE_MT", "CPUE_LBS")))
  
  
  # If 'output' is specified, return appropriate data
  if(!is.null(output)){
    if(output == "cpue"){
      return(cpue = cpue_out)
    }
    
    if(output == "bioabund"){
      # remove station with NA (2000 E-11 BB males)
      station_cpue <- station_cpue %>%
        mutate(CPUE = log(CPUE+1)) %>%
        filter(!is.na(COUNT))
      
      return(list(group_cols = groups_out,
                  cpue = station_cpue))
    }
  }
  
  # If 'output' isn't specified, default output is 'cpue'
  if(is.null(output)){
    return(cpue = cpue_out)
  }
  
}



log_bioabund <- function(crab_data = NULL,
                          species = NULL,
                          region = c("EBS", "NBS")[1],
                          district = c("ALL", "BB", "NORTH", "NS", "PRIB", "STMATT",
                                       "UNSTRAT", "E166", "W166", "166TO173")[1],
                          years = NULL,
                          sex = NULL,
                          size_min = NULL,
                          size_max = NULL,
                          crab_category = NULL,
                          female_maturity = c("morphometric", "cutline")[1],
                          shell_condition = NULL,
                          egg_condition = NULL,
                          clutch_size = NULL,
                          bin_1mm = FALSE,
                          spatial_level = c("stratum", "district", "region")[2],
                          # rm_corners = FALSE,
                          replace_retow = TRUE){
  
  
  ## ERROR: must specify just one output if length(output > 1)
  if(length(spatial_level) > 1){
    stop("Argument `spatial_level` must be of length = 1.")
  }
  
  
  # call calc_cpue()
  cpue <- log_cpue(crab_data = crab_data,
                              species = species,
                              region = region,
                              district = district,
                              years = years,
                              sex = sex,
                              size_min = size_min,
                              size_max = size_max,
                              crab_category = crab_category,
                              female_maturity = female_maturity,
                              shell_condition = shell_condition,
                              egg_condition = egg_condition,
                              clutch_size = clutch_size,
                              bin_1mm = bin_1mm,
                              replace_retow = replace_retow,
                              # rm_corners = rm_corners,
                              output = c("bioabund"))
  
  station_cpue <- cpue$cpue
  group_cols <- cpue$group_cols
  
  # make a vector of years to re-expand zero-catch years by
  year_vec <- unique(station_cpue$YEAR)
  
  
  # Identify zero-catch stations for Northern District RKC and Hair Crab, and BKC Unstratified
  if(TRUE %in% (c("NORTH", "UNSTRAT") %in% unique(station_cpue$DISTRICT))){
    # Filter out zero-catch stations if District is NORTH or UNSTRAT
    zero_catch <- station_cpue %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'REGION', 'DISTRICT', 'STRATUM', 'TOTAL_AREA', 'STATION_ID')))) %>%
      dplyr::reframe(TOTAL_CPUE = sum(CPUE)) %>%
      dplyr::mutate(REMOVE = dplyr::case_when((DISTRICT %in% c("NORTH", "UNSTRAT") & TOTAL_CPUE == 0) ~ "remove",
                                              TRUE ~ "keep"))
    
    station_cpue <- station_cpue %>%
      dplyr::left_join(., zero_catch,
                       by = c('YEAR', 'STATION_ID', 'REGION',
                              'DISTRICT', 'STRATUM', 'TOTAL_AREA')) %>%
      dplyr::filter(REMOVE == "keep") %>%
      dplyr::select(-c('REMOVE', 'TOTAL_CPUE'))
  }
  
  
  
  ## Calculate abundance and biomass
  # Sum across haul, scale abundance, biomass, and variance to strata, then sum across strata and calc CIs
  bio_abund_stratum <- station_cpue %>%
    # Scale to abundance by strata
    dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'REGION', 'DISTRICT', 'STRATUM', 'TOTAL_AREA', group_cols)))) %>%
    dplyr::reframe(MEAN_CPUE = mean(CPUE),
                   N_CPUE = dplyr::n(),
                   VAR_CPUE = (stats::var(CPUE)*(TOTAL_AREA^2))/N_CPUE,
                   SD_CPUE = sqrt(VAR_CPUE),
                   MEAN_CPUE_MT = mean(CPUE_MT),
                   N_CPUE_MT = dplyr::n(),
                   VAR_CPUE_MT = (stats::var(CPUE_MT)*(TOTAL_AREA^2))/N_CPUE_MT,
                   SD_CPUE_MT = sqrt(VAR_CPUE_MT),
                   MEAN_CPUE_LBS = mean(CPUE_LBS),
                   N_CPUE_LBS = dplyr::n(),
                   VAR_CPUE_LBS = (stats::var(CPUE_LBS)*(TOTAL_AREA^2))/N_CPUE_LBS,
                   SD_CPUE_LBS = sqrt(VAR_CPUE_LBS),
                   ABUNDANCE = (MEAN_CPUE * TOTAL_AREA),
                   ABUNDANCE_CV = (SD_CPUE/ABUNDANCE),
                   ABUNDANCE_CI = 1.96*(SD_CPUE),
                   BIOMASS_MT = (MEAN_CPUE_MT * TOTAL_AREA),
                   BIOMASS_MT_CV = (SD_CPUE_MT/BIOMASS_MT),
                   BIOMASS_MT_CI = (1.96*SD_CPUE_MT),
                   BIOMASS_LBS = (MEAN_CPUE_LBS * TOTAL_AREA),
                   BIOMASS_LBS_CV = (SD_CPUE_LBS/BIOMASS_LBS),
                   BIOMASS_LBS_CI = 1.96*(SD_CPUE_LBS),
                   N_STATIONS = length(unique(STATION_ID))) %>%
    dplyr::distinct()
  
  
  # Re-expand by year to add back in 0-catch years for Northern District RKC and Hair Crab, and BKC Unstratified
  if(TRUE %in% (c("NORTH", "UNSTRAT") %in% unique(station_cpue$DISTRICT))){
    bio_abund_stratum <- bio_abund_stratum %>%
      dplyr::right_join(tidyr::expand_grid(SEX_TEXT = sex_combos,
                                           CATEGORY = category_combos,
                                           SHELL_TEXT = shell_combos,
                                           EGG_CONDITION_TEXT = egg_combos,
                                           CLUTCH_TEXT = clutch_combos,
                                           SIZE_1MM = bin_combos,
                                           YEAR = year_vec,
                                           REGION = unique(station_cpue$REGION),
                                           DISTRICT = unique(station_cpue$DISTRICT)),
                        by = c('YEAR', 'REGION', 'DISTRICT', group_cols),
                        relationship = "many-to-many") %>%
      dplyr::select(dplyr::all_of(c("YEAR", "REGION", "DISTRICT", "STRATUM",
                                    "TOTAL_AREA", "N_STATIONS", group_cols,
                                    "MEAN_CPUE", "N_CPUE", "VAR_CPUE", "SD_CPUE",
                                    "MEAN_CPUE_MT", "N_CPUE_MT", "VAR_CPUE_MT", "SD_CPUE_MT",
                                    "MEAN_CPUE_LBS", "N_CPUE_LBS", "VAR_CPUE_LBS", "SD_CPUE_LBS",
                                    "ABUNDANCE", "ABUNDANCE_CV", "ABUNDANCE_CI",
                                    "BIOMASS_MT", "BIOMASS_MT_CV", "BIOMASS_MT_CI",
                                    "BIOMASS_LBS", "BIOMASS_LBS_CV", "BIOMASS_LBS_CI"))) %>%
      distinct()
  }
  
  if(spatial_level == "stratum"){
    # Format output dataframe
    stratum_out <- bio_abund_stratum %>%
      dplyr::mutate(SPECIES = species) %>%
      dplyr::select(dplyr::all_of(c('SPECIES', 'YEAR', 'REGION', 'DISTRICT',
                                    'STRATUM', 'TOTAL_AREA', group_cols,
                                    'ABUNDANCE', 'ABUNDANCE_CV', 'ABUNDANCE_CI',
                                    'BIOMASS_MT', 'BIOMASS_MT_CV', 'BIOMASS_MT_CI',
                                    'BIOMASS_LBS', 'BIOMASS_LBS_CV', 'BIOMASS_LBS_CI'))) %>%
      replace(is.na(.), 0)
    
    return(stratum_out)
  }
  
  
  bio_abund_district <- bio_abund_stratum %>%
    # Sum across strata
    dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'REGION', 'DISTRICT', group_cols)))) %>%
    dplyr::reframe(TOTAL_AREA = sum(TOTAL_AREA),
                   MEAN_CPUE = sum(MEAN_CPUE),
                   VAR_CPUE = sum(VAR_CPUE),
                   SD_CPUE = sqrt(VAR_CPUE),
                   N_CPUE = sum(N_CPUE),
                   MEAN_CPUE_MT = sum(MEAN_CPUE_MT),
                   VAR_CPUE_MT = sum(VAR_CPUE_MT),
                   SD_CPUE_MT = sqrt(VAR_CPUE_MT),
                   N_CPUE_MT = sum(N_CPUE_MT),
                   MEAN_CPUE_LBS = sum(MEAN_CPUE_LBS),
                   VAR_CPUE_LBS = sum(VAR_CPUE_LBS),
                   SD_CPUE_LBS = sqrt(VAR_CPUE_LBS),
                   N_CPUE_LBS = sum(N_CPUE_LBS),
                   ABUNDANCE = sum(ABUNDANCE),
                   ABUNDANCE_CV = (SD_CPUE/ABUNDANCE),
                   ABUNDANCE_CI = 1.96*(SD_CPUE),
                   BIOMASS_MT = sum(BIOMASS_MT),
                   BIOMASS_MT_CV = (SD_CPUE_MT/BIOMASS_MT),
                   BIOMASS_MT_CI = (1.96*SD_CPUE_MT),
                   BIOMASS_LBS = sum(BIOMASS_LBS),
                   BIOMASS_LBS_CV = (SD_CPUE_LBS/BIOMASS_LBS),
                   BIOMASS_LBS_CI = 1.96*(SD_CPUE_LBS),
                   N_STATIONS = sum(N_STATIONS)) %>%
    dplyr::mutate(N_STATIONS = ifelse((YEAR == 2000 & DISTRICT == "BB"), 135, N_STATIONS)) %>%
    dplyr::ungroup()
  
  if(spatial_level == "district"){
    # Format output dataframe
    district_out <- bio_abund_district %>%
      dplyr::mutate(SPECIES = species) %>%
      dplyr::select(dplyr::all_of(c('SPECIES', 'YEAR', 'REGION', 'DISTRICT', 'TOTAL_AREA', group_cols,
                                    'ABUNDANCE', 'ABUNDANCE_CV', 'ABUNDANCE_CI',
                                    'BIOMASS_MT', 'BIOMASS_MT_CV', 'BIOMASS_MT_CI',
                                    'BIOMASS_LBS', 'BIOMASS_LBS_CV', 'BIOMASS_LBS_CI'))) %>%
      replace(is.na(.), 0)
    return(district_out)
  }
  
  
  bio_abund_region <- bio_abund_district %>%
    # Sum across districts
    dplyr::group_by(dplyr::across(dplyr::all_of(c('YEAR', 'REGION', group_cols)))) %>%
    dplyr::reframe(TOTAL_AREA = sum(TOTAL_AREA),
                   MEAN_CPUE = sum(MEAN_CPUE),
                   VAR_CPUE = sum(VAR_CPUE),
                   SD_CPUE = sqrt(VAR_CPUE),
                   N_CPUE = sum(N_CPUE),
                   MEAN_CPUE_MT = sum(MEAN_CPUE_MT),
                   VAR_CPUE_MT = sum(VAR_CPUE_MT),
                   SD_CPUE_MT = sqrt(VAR_CPUE_MT),
                   N_CPUE_MT = sum(N_CPUE_MT),
                   MEAN_CPUE_LBS = sum(MEAN_CPUE_LBS),
                   VAR_CPUE_LBS = sum(VAR_CPUE_LBS),
                   SD_CPUE_LBS = sqrt(VAR_CPUE_LBS),
                   N_CPUE_LBS = sum(N_CPUE_LBS),
                   ABUNDANCE = sum(ABUNDANCE),
                   ABUNDANCE_CV = (SD_CPUE/ABUNDANCE),
                   ABUNDANCE_CI = 1.96*(SD_CPUE),
                   BIOMASS_MT = sum(BIOMASS_MT),
                   BIOMASS_MT_CV = (SD_CPUE_MT/BIOMASS_MT),
                   BIOMASS_MT_CI = (1.96*SD_CPUE_MT),
                   BIOMASS_LBS = sum(BIOMASS_LBS),
                   BIOMASS_LBS_CV = (SD_CPUE_LBS/BIOMASS_LBS),
                   BIOMASS_LBS_CI = 1.96*(SD_CPUE_LBS),
                   N_STATIONS = sum(N_STATIONS)) %>%
    dplyr::ungroup()
  
  if(spatial_level == "region"){
    # Format output dataframe
    region_out <- bio_abund_region %>%
      dplyr::mutate(SPECIES = species) %>%
      dplyr::select(dplyr::all_of(c('SPECIES', 'YEAR', 'REGION', 'TOTAL_AREA', group_cols,
                                    'ABUNDANCE', 'ABUNDANCE_CV', 'ABUNDANCE_CI',
                                    'BIOMASS_MT', 'BIOMASS_MT_CV', 'BIOMASS_MT_CI',
                                    'BIOMASS_LBS', 'BIOMASS_LBS_CV', 'BIOMASS_LBS_CI'))) %>%
      replace(is.na(.), 0)
    return(region_out)
  }
  
}
