## Purpose: To calculate (1) mean bottom temperature, (2) cold pool extent, 
##          (3) mean surface temperature, (4) 
##          from the EBS Bottom Trawl Survey time series. 
##
## NOTES:
## - using 1988+, but think about imputing temperature for missing years


## Read in setup
source("./scripts/setup.R")


## Extract haul data
haul <- tanner$haul


## Compute mean summer bottom temperature
mean_bt <- haul %>%
           filter(YEAR %in% years,
                  !HAUL_TYPE == 17) %>%
           distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
           group_by(YEAR) %>%
           summarise(summer_bt = mean(GEAR_TEMPERATURE, na.rm = T)) %>%
           right_join(expand.grid(YEAR = years))

# Plot
mean_bt %>%
  ggplot(aes(x = as.numeric(YEAR), y = summer_bt)) +
  geom_point() +
  geom_line() +
  labs(y = "Bottom Temperature (C)", x = "Year") +
  geom_hline(aes(yintercept = mean(summer_bt, na.rm = TRUE)), linetype = 5) +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave("./figures/bottom_temp.png")



## Compute cold pool areal extent
cp_extent <- haul %>%
             filter(YEAR %in% years,
                    !HAUL_TYPE == 17,
                    !(STATION_ID %in% corners)) %>%
             distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
             group_by(YEAR) %>%
             summarise(cp_extent = sum(GEAR_TEMPERATURE < 2, na.rm = T) * 401) %>%
             right_join(expand.grid(YEAR = years))

# Plot
cp_extent %>%
  ggplot(aes(x = as.numeric(YEAR), y = cp_extent)) +
  geom_point() +
  geom_line() +
  labs(y = "Cold Pool Extent (nmi2)", x = "Year") +
  geom_hline(aes(yintercept = mean(cp_extent, na.rm = TRUE)), linetype = 5) +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave("./figures/coldpool_extent.png")



## Compute mean summer surface temperature
mean_st <- haul %>%
           filter(YEAR %in% years,
                  !HAUL_TYPE == 17) %>%
           distinct(YEAR, STATION_ID, SURFACE_TEMPERATURE) %>%
           group_by(YEAR) %>%
           summarise(summer_st = mean(SURFACE_TEMPERATURE, na.rm = T)) %>%
           right_join(expand.grid(YEAR = years))

# Plot
mean_st %>%
  ggplot(aes(x = as.numeric(YEAR), y = summer_st)) +
  geom_point() +
  geom_line() +
  labs(y = "Surface Temperature (C)", x = "Year") +
  geom_hline(aes(yintercept = mean(summer_st, na.rm = TRUE)), linetype = 5) +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave("./figures/surface_temp.png")





## Combine indices and save output
full_join(mean_bt, cp_extent, mean_st) %>%
  write_csv("./outputs/temp_coldpool.csv")








