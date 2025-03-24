## Purpose: To calculate mean bottom temperature and cold pool extent indices
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
           filter(YEAR >= 1988,
                  HAUL_TYPE == 3) %>%
           distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
           group_by(YEAR) %>%
           summarise(summer_bt = mean(GEAR_TEMPERATURE, na.rm = T)) %>%
           right_join(expand.grid(YEAR = c(1988:current_year)))

# Plot
mean_bt %>%
  ggplot(aes(x = as.numeric(YEAR), y = summer_bt)) +
  geom_point() +
  geom_line() +
  labs(y = "Bottom Temperature (C)", x = "Year") +
  geom_hline(aes(yintercept = mean(summer_bt, na.rm = TRUE)), linetype = 5) +
  xlim(1988, 2024) +
  theme_bw()
ggsave("./figures/bottom_temp.png")



## Compute cold pool areal extent
cp_extent <- haul %>%
             filter(YEAR >= 1988,
                    HAUL_TYPE == 3,
                    !(STATION_ID %in% corners)) %>%
             distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
             group_by(YEAR) %>%
             summarise(cp_extent = sum(GEAR_TEMPERATURE < 2, na.rm = T) * 401) %>%
             right_join(expand.grid(YEAR = c(1988:current_year)))

# Plot
cp_extent %>%
  ggplot(aes(x = as.numeric(YEAR), y = cp_extent)) +
  geom_point() +
  geom_line() +
  labs(y = "Cold Pool Extent (nmi2)", x = "Year") +
  geom_hline(aes(yintercept = mean(cp_extent, na.rm = TRUE)), linetype = 5) +
  xlim(1988, 2024) +
  theme_bw()
ggsave("./figures/coldpool_extent.png")



## Combine indices and save output
full_join(mean_bt, cp_extent) %>%
  write_csv("./outputs/temp_coldpool.csv")








