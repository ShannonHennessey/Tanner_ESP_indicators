## Purpose: To summarize benthic predator and Pacific cod mean CPUE across years 
##          in core tanner habitat
##
## NOTES:
## - revise these to a spatial overlap metric?
## - predator biomass vs. abundance??
## - SIZE CUTOFF - crab and preds (maturity probably doesn't matter....)
## - try to refine suite of predators....


## Load packages
library(crabpack)
library(tidyverse)
# library(mgcv)
library(dispRity)
library(rnaturalearth)
library(ggnewscale)
library(stars)
library(sf)
library(terra)
library(ENMTools)


## Read in setup
source("./scripts/setup.R")


## Function to calculate Hellinger's distance
norm_vec <- function(x) sqrt(sum(x^2))

hell_dist <- function (p, q, from, to, n = 1024) {
  P <- density(p, kernel = "gaussian", from = from, to = to, n = n)
  p <- P$y
  p <- p / sum(p)
  Q <- density(q, kernel = "gaussian", from = from, to = to, n = n)
  q <- Q$y
  q <- q / sum(q)
  hd <- norm_vec(sqrt(p) - sqrt(q)) / sqrt(2)
  hd
}


## Read in data 
# Load groundfish data queried from Oracle -- from "get_gf_data.R" script
pred <- read.csv(paste0(data_dir, "gf_cpue_timeseries.csv")) %>%
        rename(STATION_ID = STATION)

# Load tanner core area stations -- from "get_core_area.R" script
core_stations <- read_csv("./outputs/tanner_area_cpue50.csv", show_col_types = FALSE) %>% 
                 pull(STATION_ID)


## Calculate mean CPUE for each predator guild across years 
# Specifying each species here because stomach contents/diets were validated for most
# and included if assumed to be benthic predator on crab juv/adult stages 

# Define guilds by species code
sab_hal <- c(20510, 10120)
pcod <- c(21720, 21722)
skate <- c(420,435,440,455,471,472,480,460,485)
flatfish <- c(10220,10115,10130,10140,10120,10260,10261,10210)
sculpin <- c(21347,21348,21368,21370,21388,21420,21311,21315,21390,21438,21371)
eelpout <- c(24184, 24191, 24185)
wolffish <- c(20320, 20322)
octopus <- c(78010, 78012, 78403)

flathead_sole <- c(10130)
rock_sole <- c(10260, 10261)
yellowfin_sole <- c(10210)


guilds <- c(#"sab_hal",  "skate", "flatfish", "sculpin", "eelpout", "wolffish", "octopus",
            "pcod", "flathead_sole", "rock_sole", "yellowfin_sole")

stations <- pred %>%
            select(YEAR, STATION_ID) %>%
            filter(STATION_ID %in% core_stations, 
                   YEAR %in% years) %>%
            distinct()

# Calculate mean CPUE by guild and year  
ben_pred <- pred %>%
            mutate(GUILD = case_when(# SPECIES_CODE %in% sab_hal ~ "sab_hal",
                                     SPECIES_CODE %in% pcod ~ "pcod",
                                     # SPECIES_CODE %in% skate ~ "skate",
                                     # SPECIES_CODE %in% flatfish ~ "flatfish",
                                     # SPECIES_CODE %in% sculpin ~ "sculpin",
                                     # SPECIES_CODE %in% eelpout ~ "eelpout",
                                     # SPECIES_CODE %in% wolffish ~ "wolffish",
                                     # SPECIES_CODE %in% octopus ~ "octopus",
                                     SPECIES_CODE %in% flathead_sole ~ "flathead_sole",
                                     SPECIES_CODE %in% rock_sole ~ "rock_sole",
                                     SPECIES_CODE %in% yellowfin_sole ~ "yellowfin_sole",
                                     TRUE ~ NA)) %>%             
            filter(STATION_ID %in% core_stations, 
                   YEAR %in% years,
                   !is.na(GUILD)) %>%   
            # station-level cpue by guild
            group_by(YEAR, STATION_ID, GUILD) %>%
            summarise(CPUE_KGKM2 = sum(CPUE_KGKM2)) %>%
            # add in 0-catch stations by guild
            right_join(., expand_grid(stations, GUILD = guilds)) %>% 
            arrange(YEAR, STATION_ID, GUILD) %>%
            mutate(CPUE_KGKM2 = replace_na(CPUE_KGKM2, 0)) %>%
            # annual mean cpue by guild
            group_by(YEAR, GUILD) %>%
            summarise(CPUE_KGKM2 = mean(CPUE_KGKM2)) %>%
            right_join(., expand.grid(YEAR = years, 
                                      GUILD = guilds)) %>%
            arrange(YEAR, GUILD) %>%
            mutate(CPUE_KGKM2 = ifelse(YEAR == 2020, NA, CPUE_KGKM2))

# Plot 
guild_plot <- ben_pred %>%
              # pivot_longer(c(2:9), names_to = "pred_guild", values_to = "CPUE_KGKM2") %>%
              ggplot(aes(x = YEAR, y = CPUE_KGKM2, group = factor(GUILD))) +
              geom_point(aes(colour = GUILD)) +
              geom_line(aes(colour = GUILD)) +
              labs(y = "Benthic Predator CPUE (kg/km2)", x = "Year") +
              theme_bw() +
              theme(legend.title = element_blank())
guild_plot # biomass is really dominated by YFS

guild_facet <- ben_pred %>%
               # pivot_longer(c(2:9), names_to = "pred_guild", values_to = "CPUE_KGKM2") %>%
               ggplot(aes(x = YEAR, y = CPUE_KGKM2)) +
               geom_point() +
               geom_line() +
               labs(y = "CPUE (kg/km2)", x = "Year") +
               theme_bw() +
               theme(legend.title = element_blank()) +
               facet_wrap(~GUILD, scales = "free_y", ncol = 8)
ggsave(paste0(fig_dir, "benthic_predator_facet.png"), guild_facet,
       height = 3, width = 15)

ben_pred %>%
  group_by(YEAR) %>%
  mutate(total_pred = sum(CPUE_KGKM2)) %>%
  ggplot(aes(x = YEAR, y = total_pred)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(total_pred, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(total_pred, na.rm = TRUE) - sd(total_pred, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(total_pred, na.rm = TRUE) + sd(total_pred, na.rm = TRUE)), color = "green4") +
  labs(y = "Benthic Predator\nCPUE (kg/km2)", x = "Year") +
  theme_bw() +
  theme(legend.title = element_blank()) 
ggsave(paste0(fig_dir, "benthic_predator_density.png"),
       height = 2, width = 6)

# Plot Pacific cod
pcod_plot <- ben_pred %>%
             ggplot(aes(x = YEAR, y = pcod)) +
             geom_point() +
             geom_line() +
             labs(y = "Pacific Cod CPUE (kg/km2)", x = "Year") +
             theme_bw()
ggsave(paste0(fig_dir, "pcod_density.png"), pcod_plot,
       height = 6, width = 10)


# Write .csv output of benthic predator density
ben_pred %>%
  group_by(YEAR) %>%
  summarise(total_pred = sum(CPUE_KGKM2)) %>%
  rename(year = YEAR) %>%
  write.csv("./outputs/benthic_predator_density.csv", row.names = FALSE)




## OVERLAP METRICS -------------------------------------------------------------

## maybe just plot yearly maps for each species overlap with tanner to look at this visually?
## help give some sort of expectation as to how these 'overlaps' are supposed to look...

## calculate CPUE for benthic predators just at station level
# define stations sampled each year
stations <- pred %>%
            select(YEAR, STATION_ID) %>%
            filter(#STATION_ID %in% core_stations, 
                   YEAR %in% years) %>%
            distinct() %>%
            filter(!STATION_ID %in% corners) # remove corner stations



# set metric (so can change easily in pipe below)
metric <- "CPUE_KGKM2"

ben_pred <- pred %>%
            mutate(GUILD = case_when(SPECIES_CODE %in% sab_hal ~ "sab_hal",
                                     SPECIES_CODE %in% pcod ~ "pcod",
                                     SPECIES_CODE %in% skate ~ "skate",
                                     SPECIES_CODE %in% flatfish ~ "flatfish",
                                     SPECIES_CODE %in% sculpin ~ "sculpin",
                                     SPECIES_CODE %in% eelpout ~ "eelpout",
                                     SPECIES_CODE %in% wolffish ~ "wolffish",
                                     SPECIES_CODE %in% octopus ~ "octopus",
                                     TRUE ~ NA)) %>%             
            filter(#STATION_ID %in% core_stations, 
                   STATION_ID %in% stations$STATION_ID,
                   YEAR %in% years,
                   !is.na(GUILD)) %>%   
            # station-level cpue by guild
            group_by(YEAR, STATION_ID, GUILD) %>%
            summarise(COUNT = sum(COUNT),
                      CPUE_KGKM2 = sum(CPUE_KGKM2),
                      CPUE_NOKM2 = sum(CPUE_NOKM2)) %>%
            # add in 0-catch stations by guild
            right_join(., expand_grid(stations, GUILD = guilds)) %>% 
            arrange(YEAR, STATION_ID, GUILD) %>%
            mutate(COUNT = replace_na(COUNT, 0),
                   CPUE_KGKM2 = replace_na(CPUE_KGKM2, 0),
                   CPUE_NOKM2 = replace_na(CPUE_NOKM2, 0))



## Pull size at 50% probability of terminal molt
# Assign static mean cutline to missing years:
# 103.5mm population, 110mm E166, 99mm W166
mat_size <- get_male_maturity(species = "TANNER", 
                              region = "EBS")$model_parameters %>% 
            select(-c("A_EST", "A_SE")) %>%
            rename(MAT_SIZE = B_EST, 
                   STD_ERR = B_SE) %>%
            right_join(., expand_grid(YEAR = years,
                                      SPECIES = "TANNER", 
                                      REGION = "EBS",
                                      DISTRICT = c("ALL", "E166", "W166"))) %>%
            mutate(MAT_SIZE = case_when(DISTRICT == "ALL" & is.na(MAT_SIZE) ~ 103.5, 
                                        DISTRICT == "E166" & is.na(MAT_SIZE) ~ 110, 
                                        DISTRICT == "W166" & is.na(MAT_SIZE) ~ 99, 
                                        TRUE ~ MAT_SIZE))


## Compute station-level CPUE by size-sex category
# Assign maturity to specimen data; calculate CPUE
cpue <- tanner$specimen %>% 
        left_join(., mat_size) %>%
        mutate(CATEGORY = case_when((SEX == 1 & SIZE >= MAT_SIZE) ~ "mature_male",
                                    (SEX == 1 & SIZE < MAT_SIZE) ~ "immature_male",
                                    (SEX == 2 & CLUTCH_SIZE >= 1) ~ "mature_female",
                                    (SEX == 2 & CLUTCH_SIZE == 0) ~ "immature_female",
                                    TRUE ~ NA)) %>%
        filter(YEAR %in% years,
               !is.na(CATEGORY)) %>%
        group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE, AREA_SWEPT, CATEGORY) %>%
        summarise(COUNT = sum(SAMPLING_FACTOR)) %>%
        pivot_wider(names_from = CATEGORY, values_from = COUNT) %>%
        mutate(population = sum(immature_male, mature_male, immature_female, mature_female, na.rm = T)) %>%
        pivot_longer(c(6:10), names_to = "CATEGORY", values_to = "COUNT") %>%
        filter(CATEGORY != "NA") %>%
        mutate(COUNT = replace_na(COUNT, 0),
               CPUE = COUNT / AREA_SWEPT) %>%
        ungroup() %>%
        select(-COUNT) %>%
        pivot_wider(names_from = CATEGORY, values_from = CPUE) %>%
        right_join(., tanner$haul %>% 
                        select(YEAR, STATION_ID) %>% 
                        filter(YEAR %in% years)) %>%
        mutate(immature_male = replace_na(immature_male, 0),
               mature_male = replace_na(mature_male, 0),
               immature_female = replace_na(immature_female, 0),
               mature_female = replace_na(mature_female, 0),
               population = replace_na(population, 0),
               mature = mature_male + mature_female,
               immature = immature_male + immature_female) %>%
        filter(STATION_ID %in% stations$STATION_ID) # remove corner stations
        

## Bind pred/tanner cpue
cpue_all <- cpue %>%
            select(-LATITUDE, -LONGITUDE, -AREA_SWEPT) %>%
            left_join(., ben_pred) %>%
            left_join(., tanner$haul %>% 
                          select(YEAR, STATION_ID, MID_LATITUDE, MID_LONGITUDE) %>%
                          rename(LATITUDE = MID_LATITUDE,
                                 LONGITUDE = MID_LONGITUDE)) 





## loop over crab categories to calculate overlaps
crab_categories <- c("population", "mature", "immature")#, "immature_male", "mature_male", "immature_female", "mature_female")

# Set crs, years, metric (so can change easily in pipe below)
map.crs <- "EPSG:3338"
years <- c(1988:2019, 2021:2024)
metric <- "CPUE_KGKM2"

corr_df <- c()

for(c in 1:length(crab_categories)){
  # set crab category
  crab <- crab_categories[c]
  
  # calculate overlap coefficient by year
  for(i in 1:length(years)){

    # filter and rescale data
    dat <- cpue_all %>%
           filter(YEAR == years[i]) %>%
           select(all_of(c("YEAR", "STATION_ID", "LATITUDE", "LONGITUDE", 
                           crab, "GUILD", metric))) %>%
           group_by(GUILD, YEAR) %>%
           # rescale for warren's I and Schoener's D
           mutate("{crab}" := scales::rescale(get(crab)), 
                  "{metric}" := scales::rescale(get(metric))) %>%
           ungroup() 
    
    # loop over guilds
    for(g in 1:length(guilds)){
      guild <- guilds[g]  
         
      crab_rast <- dat %>%
                   filter(GUILD == guild) %>%
                   rename(value = crab) %>%
                   select(LATITUDE, LONGITUDE, value) %>%
                   sf::st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), 
                                crs = map.crs) %>%
                   stars::st_rasterize(.) %>%
                   terra::rast(.) # SpatRaster
                   # raster::raster(.) # RasterLayer  
       
      guild_rast <- dat %>%
                    filter(GUILD == guild) %>%
                    rename(value = metric) %>%
                    select(LATITUDE, LONGITUDE, value) %>%
                    sf::st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), 
                                 crs = map.crs) %>%
                    stars::st_rasterize(.) %>%
                    terra::rast(.) # SpatRaster
                    # raster::raster(.) # RasterLayer  
      
      # skip year if guild data == 0
      # skates, 1989:1995, 1997
      if(years[i] %in% c(1989:1995, 1997) & guild == "skate"){
        out <- data.frame(year = years[i],
                          crab = crab,
                          guild = guild,
                          corr = NA,
                          I = NA,
                          SD = NA,
                          HD = NA)
        
        corr_df <- rbind(corr_df, out)
        next
      }
      
      
      # calculate overlap metrics
      overlap <- ENMTools::raster.overlap(crab_rast, guild_rast)
      
      
      # set bins for Hellinger's distance
      n_bins <- 1024 #arbitrary
      bins <- seq(floor(min(na.omit(c(values(crab_rast), values(guild_rast))))), 
                  ceiling(max(c(na.omit(values(crab_rast), values(guild_rast))))), length = n_bins)
      
      # calculate Hellinger's distance
      hd <- hell_dist(scale(na.omit(values(crab_rast))), scale(na.omit(values(guild_rast))), 
                      from = bins[1], to = bins[n_bins])
      
      # compile metrics
      out <- data.frame(year = years[i],
                       crab = crab_categories[c],
                       guild = guilds[g],
                       corr = round(overlap$rank.cor, 2),
                       I = round(overlap$I, 2),
                       SD = round(overlap$D, 2),
                       HD = hd)
      
      corr_df <- rbind(corr_df, out)
    } # end guild loop
  } # end year loop
  
  # overlap_coeff <- cpue_all %>%
  #                  filter(!YEAR == 2020) %>%
  #                  select(all_of(c("YEAR", "STATION_ID", crab, "GUILD", metric))) %>%
  #                  group_by(GUILD, YEAR) %>%
  #                  mutate("{crab}" := BBmisc::normalize(get(crab), range = c(0, 1)), 
  #                         "{metric}" := BBmisc::normalize(get(metric), range = c(0, 1))) %>%
  #                  summarise(COEFF = dispRity::bhatt.coeff(get(crab), get(metric))) %>%
  #                  right_join(., expand_grid(YEAR = c(1988:2024),
  #                                            GUILD = guilds)) %>%
  #                  arrange(YEAR) %>%
  #                  mutate(MEAN = mean(COEFF, na.rm = TRUE))
  # 
  # 
  # ## Plot coefficients
  # coeff_plot <- ggplot(overlap_coeff, aes(x = YEAR, y = COEFF, group = factor(GUILD))) +
  #               geom_point(aes(colour = GUILD)) +
  #               geom_line(aes(colour = GUILD)) +
  #               geom_hline(aes(yintercept = MEAN), linetype = 5) +
  #               labs(y = "Overlap coefficient", x = "Year") +
  #               theme_bw() +
  #               theme(legend.title = element_blank()) +
  #               facet_wrap(~GUILD, scales = "free")
  # # coeff_plot
  # ggsave(paste0(fig_dir, "exploratory/predator_overlap_", crab, ".png"), 
  #        coeff_plot, height = 6, width = 10)
  # 
  # 
  # 
  #   
  #   ## maps 
  #   guild_overlap <- ggplot(ne_countries(scale = "medium", returnclass = "sf")) +
  #                    geom_sf() +
  #                    geom_point(data = cpue_all %>% filter(GUILD == guild), 
  #                               aes(x = LONGITUDE, y = LATITUDE, 
  #                                   size = get(crab), 
  #                                   shape = get(crab) == 0,
  #                                   color = get(crab) == 0), 
  #                               alpha = 1) +                  
  #                    scale_color_manual(values = c('TRUE' = "gray", 'FALSE' = "darkblue"), guide = "none") +
  #                    scale_shape_manual(values = c('TRUE' = 4, 'FALSE' = 19), guide = "none") +        
  #                    new_scale("color") +
  #                    geom_point(data = cpue_all %>% filter(GUILD == guild),
  #                               aes(x = LONGITUDE, y = LATITUDE,
  #                                   size = get(metric),
  #                                   shape = get(metric) == 0,
  #                                   color = get(metric) == 0), 
  #                               alpha = 0.7) +
  #                    scale_color_manual(values = c('TRUE' = "gray", 'FALSE' = "red"), guide = "none") +
  #                    coord_sf(xlim = c(-180, -157), ylim = c(54, 62.5), expand = FALSE) +
  #                    labs(title = ) +
  #                    theme_bw() +
  #                    theme(legend.position = "none") +
  #                    geom_text(data = overlap_coeff %>% filter(GUILD == guild),
  #                              mapping = aes(x = -Inf, y = -Inf, label = round(COEFF, digits = 2)),
  #                              hjust = -0.1,
  #                              vjust = -1) +
  #                    facet_wrap(~YEAR)
  #   ggsave(paste0(fig_dir, "exploratory/", guild, "_overlap_", crab, ".png"), guild_overlap,
  #          height = 20, width = 30)
  # } # end guild loop
} # end crab category loop


## plot overlap metric timeseries by crab category and guild
# corr, I, SD, HD
correlations <- corr_df %>%
                group_by(crab, guild) %>%
                mutate(MEAN = mean(corr, na.rm = TRUE)) %>%
                right_join(., expand_grid(year = c(1988:2024),
                                          crab = crab_categories, 
                                          guild = guilds))

## do some sort of test to see if maturity scores sig diff from population?
                
ggplot(correlations, aes(x = year, y = (HD), color = as.factor(crab))) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = (MEAN), color = as.factor(crab)), linetype = 5) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~crab + guild, scales = 'free', ncol = 8)
ggsave(paste0(fig_dir, "exploratory/pred_overlap_HD.png"), height = 6, width = 15)




# ## Write .csv for stations in 50th percentile of avg CPUE  
# write.csv(cpue50_core %>% 
#             filter(MEAN_CPUE > quantile(MEAN_CPUE, 0.50)) %>%
#             select(-QUANTILE), 
#           "./outputs/tanner_area_cpue50.csv", row.names = FALSE)
# 



# incorporate a lag??




