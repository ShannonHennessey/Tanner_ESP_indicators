## Purpose: To calculate (1) mean bottom temperature, (2) cold pool extent, 
##          (3) mean surface temperature from the EBS Bottom Trawl Survey timeseries,
##          and (4) mean winter Aleutian Low - Beaufort Sea Anticyclone 
##          (https://psl.noaa.gov/data/timeseries/ALBSA/) using 850mb height from 
##          NCEP R1 4 pts.
##
## NOTES:
## - using 1988+, but think about imputing temperature for missing years
## - bottom temperature: positive? for both juveniles and larvae maybe? increased gonadogenesis/enbryo incubation
## - surface temperature - warm = good, favorable production of copepod nauplii?
##
## THINK ABOUT: 
## - bottom temp, cold pool extent, juv temp occupancy all highly correlated...
##   - does it make sense to use all 3 in DFA, or just pick one? (prob juv temp)
##   - only thing is eg. 1999, lg cold pool but avg juv temp; 2012 lg cold pool, very low juv temp
##     - thoughts?? 1999 avoided cold pool, 2012 didn't avoid....
##     - immature males (and kind of mature males) had low area occupied in those years...but not females
##     - 2012 had much higher density of immature males (~2x)
##
## ALPI: 3-year running averages?? -- no, because acting on larval stage
##   - Nov-Jan? Newman et al. 2016 PDO revisited - peak months for ALPI association


## Read in setup
source("./scripts/setup.R")


## Extract haul data
haul <- tanner$haul
haul_years <- c(1975:current_year)

## Compute mean summer bottom temperature
mean_bt <- haul %>%
           filter(YEAR >= 1975,
                  !HAUL_TYPE == 17) %>%
           distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
           group_by(YEAR) %>%
           summarise(summer_bt = mean(GEAR_TEMPERATURE, na.rm = T)) %>%
           right_join(expand.grid(YEAR = haul_years))

# Plot
mean_bt %>%
  ggplot(aes(x = as.numeric(YEAR), y = summer_bt)) +
  geom_point() +
  geom_line() +
  labs(y = "Bottom Temperature (C)", x = "Year") +
  geom_hline(aes(yintercept = mean(summer_bt, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(summer_bt, na.rm = TRUE) - sd(summer_bt, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(summer_bt, na.rm = TRUE) + sd(summer_bt, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "bottom_temp.png"), height = 2, width = 6)




## Compute cold pool areal extent
cp_extent <- haul %>%
             filter(YEAR >= 1975,
                    !HAUL_TYPE == 17,
                    !(STATION_ID %in% corners)) %>%
             distinct(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
             group_by(YEAR) %>%
             summarise(cp_extent = sum(GEAR_TEMPERATURE < 2, na.rm = T) * 401) %>%
             right_join(expand.grid(YEAR = haul_years))

# Plot
cp_extent %>%
  ggplot(aes(x = as.numeric(YEAR), y = cp_extent)) +
  geom_point() +
  geom_line() +
  labs(y = "Cold Pool Extent (nmi2)", x = "Year") +
  geom_hline(aes(yintercept = mean(cp_extent, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(cp_extent, na.rm = TRUE) - sd(cp_extent, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(cp_extent, na.rm = TRUE) + sd(cp_extent, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "coldpool_extent.png"), height = 2, width = 6)



## Compute mean summer surface temperature
mean_st <- haul %>%
           filter(YEAR >= 1975,
                  !HAUL_TYPE == 17) %>%
           distinct(YEAR, STATION_ID, SURFACE_TEMPERATURE) %>%
           group_by(YEAR) %>%
           summarise(summer_st = mean(SURFACE_TEMPERATURE, na.rm = T)) %>%
           right_join(expand.grid(YEAR = haul_years))

# Plot
mean_st %>%
  ggplot(aes(x = as.numeric(YEAR), y = summer_st)) +
  geom_point() +
  geom_line() +
  labs(y = "Surface Temperature (C)", x = "Year") +
  geom_hline(aes(yintercept = mean(summer_st, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(summer_st, na.rm = TRUE) - sd(summer_st, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(summer_st, na.rm = TRUE) + sd(summer_st, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "surface_temp.png"), height = 2, width = 6)




## Compute mean winter Aleutian Low index
mean_AL <- read.csv(paste0(data_dir, "albsa_daily.csv")) %>% 
           filter(month %in% c(12, 1, 2, 3)) %>% 
           mutate(year = ifelse(month == 12, year + 1, year)) %>%
           filter(year >= 1975) %>%
           group_by(year) %>%
           summarize(mean_AL = mean(ALBSA)) %>%
           rename(YEAR = year)

# Plot
mean_AL %>%
  ggplot(aes(x = YEAR, y = mean_AL)) +
  geom_point() +
  geom_line() +
  labs(y = "ALPI", x = "Year") +
  geom_hline(aes(yintercept = mean(mean_AL, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(mean_AL, na.rm = TRUE) - sd(mean_AL, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(mean_AL, na.rm = TRUE) + sd(mean_AL, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "ALPI.png"), height = 2, width = 6)



## Combine indices and save output
full_join(mean_bt, cp_extent) %>%
  full_join(., mean_st) %>%
  full_join(., mean_AL) %>%
  arrange(YEAR) %>% 
  rename(year = YEAR) %>%
  write_csv("./outputs/temp_coldpool_AL.csv")








