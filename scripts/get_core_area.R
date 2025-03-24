## Purpose: To define tanner crab core area as stations with 50th percentile of CPUE
##
## NOTES:
## - include corner stations in the quantile estimation or not??
##   there might be some spatial considerations, because not on even grid...
## - should core area be a yearly metric? (ie. time-varying??)


## Read in setup
source("./scripts/setup.R")


## Calculate station-level CPUE
cpue <- crabpack::calc_cpue(crab_data = tanner,
                            species = "TANNER",
                            years = c(1988:2024))


## Identify stations in 50-100 CPUE percentile range
cpue50_core <- cpue %>%
               filter(!STATION_ID %in% corners) %>% # remove corner stations
               group_by(STATION_ID) %>%
               summarise(MEAN_CPUE = mean(CPUE)) %>%
               mutate(QUANTILE = ifelse(MEAN_CPUE > quantile(MEAN_CPUE, 0.50), 1, 0)) %>%
               # filter(MEAN_CPUE > quantile(MEAN_CPUE, 0.50)) %>% # 187 stations, 174 excluding corners
               ungroup() %>%
               left_join(., stations %>% select(STATION_ID, LATITUDE, LONGITUDE))


## Plot for visual check
ggplot(ne_countries(scale = "medium", returnclass = "sf")) +
  geom_sf() +
  geom_point(data = cpue50_core, aes(x = LONGITUDE, y = LATITUDE, 
                                     size = MEAN_CPUE/1000, color = as.factor(QUANTILE))) +
  coord_sf(xlim = c(-180, -157), ylim = c(54, 62.5), expand = FALSE) +
  theme_bw() +
  theme(legend.position = "none")
ggsave("./figures/tanner_core_area.png", height = 6, width = 10)


## Write .csv for stations in 50th percentile of avg CPUE  
cpue50_core %>% 
  filter(MEAN_CPUE > quantile(MEAN_CPUE, 0.50)) %>%
  select(-QUANTILE) %>%
  write.csv("./outputs/tanner_area_cpue50.csv", row.names = FALSE)



