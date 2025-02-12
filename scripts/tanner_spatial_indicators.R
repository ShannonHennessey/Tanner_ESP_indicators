## Purpose: To calculate (a) Tanner Crab center of abundance (lat & lon) in EBS 
##          by size/sex category, and (b) Area Occupied (D95) - area of stations 
##          that make up 95% of the cumulative Tanner cpue
##
## Author: Shannon Hennessey; adapted from Erin Fedewa's snow crab indicaator
##
## NOTES:
## - assign male maturity by actual proportion in a given size bin, not juse a cutline?
##   ...although harder to implement for years with no cutline...
## - need to figure out what to do when missing stations in a year (see notes below for D95)

## Load packages
library(crabpack)
library(tidyverse)


## Pull Tanner specimen data
tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS")

## Pull size at 50% probability of terminal molt
# Assign static mean cutline to missing years:
# 103.5mm population, 110mm E166, 99mm W166
mat_size <- get_male_maturity(species = "TANNER", 
                              region = "EBS")$model_parameters %>% 
            select(-c("A_EST", "A_SE")) %>%
            rename(MAT_SIZE = B_EST, 
                   STD_ERR = B_SE) %>%
            right_join(., expand_grid(YEAR = c(1975:2024),
                                      SPECIES = "TANNER", 
                                      REGION = "EBS",
                                      DISTRICT = c("ALL", "E166", "W166"))) %>%
            mutate(MAT_SIZE = case_when(DISTRICT == "ALL" & is.na(MAT_SIZE) ~ 103.5, 
                                        DISTRICT == "E166" & is.na(MAT_SIZE) ~ 110, 
                                        DISTRICT == "W166" & is.na(MAT_SIZE) ~ 99, 
                                        TRUE ~ MAT_SIZE))

## Assign maturity to specimen data; calculate CPUE
cpue <- tanner$specimen %>% 
        left_join(., mat_size) %>%
        mutate(CATEGORY = case_when((SEX == 1 & SIZE >= MAT_SIZE) ~ "mature_male",
                                    (SEX == 1 & SIZE < MAT_SIZE) ~ "immature_male",
                                    (SEX == 2 & CLUTCH_SIZE >= 1) ~ "mature_female",
                                    (SEX == 2 & CLUTCH_SIZE == 0) ~ "immature_female",
                                    TRUE ~ NA)) %>%
        filter(YEAR >= 1988,
               !is.na(CATEGORY)) %>%
        group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE, AREA_SWEPT, CATEGORY) %>%
        summarise(COUNT = round(sum(SAMPLING_FACTOR))) %>%
        pivot_wider(names_from = CATEGORY, values_from = COUNT) %>%
        mutate(population = sum(immature_male, mature_male, immature_female, mature_female, na.rm = T)) %>%
        pivot_longer(c(6:10), names_to = "CATEGORY", values_to = "COUNT") %>%
        filter(CATEGORY != "NA") %>%
        mutate(COUNT = replace_na(COUNT, 0),
               CPUE = COUNT / AREA_SWEPT) %>%
        ungroup() 



## Compute Tanner crab center of abundance by size/sex
# Define corner stations
corner <- c("GF1918", "GF2019", "GF2120", "GF2221",
            "HG1918", "HG2019", "HG2120", "HG2221",
            "IH1918", "IH2019", "IH2120", "IH2221",
            "JI1918", "JI2019", "JI2120", "JI2221",
            "ON2524", "ON2625",
            "PO2423", "PO2524", "PO2625", "PO2726",
            "QP2423", "QP2524", "QP2625", "QP2726")
n_haul <- tanner_haul %>% group_by(YEAR) %>% filter(!STATION_ID %in% corner) %>% summarise(N_STATION = n())


COD <- cpue %>%
       filter(!STATION_ID %in% corner) %>% # exclude corner stations
       group_by(YEAR, CATEGORY) %>%
       summarise(LAT_COD = weighted.mean(LATITUDE, w = CPUE),
                 LON_COD = weighted.mean(LONGITUDE, w = CPUE)) %>%
       right_join(., expand_grid(YEAR = c(1989:2024),
                                 CATEGORY = unique(cpue$CATEGORY))) %>%
       arrange(YEAR)


# Plot latitude centroid
lat_cod <- ggplot(data = COD %>% filter(CATEGORY != "population"),
                  aes(x = YEAR, y = LAT_COD, group = CATEGORY, color = CATEGORY)) +
           geom_point() +
           geom_line() +
           labs(x = "", y = expression(paste("Center of Abundance (", degree, "Latitude)"))) +           
           theme_bw() +
           theme(legend.title = element_blank()) 

ggsave("./figures/tanner_centroid_latitude.png", lat_cod,
       height = 6, width = 10)


# Plot longitude centroid
lon_cod <- ggplot(data = COD %>% filter(CATEGORY != "population"),
                  aes(x = YEAR, y = LON_COD, group = CATEGORY, color = CATEGORY)) +
           geom_point() +
           geom_line() +
           labs(x = "", y = expression(paste("Center of Abundance (", degree, "Longitude)"))) +           
           theme_bw() +
           theme(legend.title = element_blank()) 

ggsave("./figures/tanner_centroid_longitude.png", lon_cod,
       height = 6, width = 10)


# Write output for COD indicator     
COD %>%
  pivot_wider(names_from = "CATEGORY", values_from = c("LAT_COD", "LON_COD")) %>%
  write.csv("./outputs/tanner_centroid_abundance.csv", row.names = FALSE)




## Compute D95 by each size and sex category 
# i.e. the number of stations contributing to 95% of cumulative CPUE
# calc from standardized timeseries (1989+)
# *note that 1992 is problematic given missing station
#  -- there are several years where not all stations were sampled...
#  -- not sure how we would deal with this??

# function to compute D95
f_d95_est <- function(x){
  x %>%
    arrange(-CPUE) %>% # sort by cpue (large:small)
    mutate(prop_cpue = CPUE/sum(CPUE),  # calculate the proportion of total cpue for each station
           cum_cpue = cumsum(prop_cpue)) %>%  
    filter(cum_cpue <= 0.95) %>% # T if in d95, F if not
    count() %>%
    mutate(d95 = (n + 1) * 401) %>% # add 1 station to n to push over 95%, multiply by 401 nm
    pull(d95)
}

# Estimate d95
d95 <- cpue %>%
       filter(!STATION_ID %in% corner) %>% # exclude corner stations
       nest(data = c(-YEAR, -CATEGORY)) %>%
       mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
       unnest(cols = c(data)) %>%
       group_by(YEAR, CATEGORY) %>%
       summarise(CPUE = sum(COUNT) / sum(AREA_SWEPT), # add a column for total cpue of each group in each year
                 d95 = mean(d95)) %>% # take 'mean' just to get one value (they are all the same)
       right_join(., expand_grid(YEAR = c(1989:2024),
                                 CATEGORY = unique(cpue$CATEGORY))) %>%
       arrange(YEAR)

# Plot
d95_plot <- ggplot(data = d95 %>% filter(CATEGORY != "population"), 
                   aes(x = YEAR, y = d95, group = CATEGORY, color = CATEGORY)) +
            geom_point() +
            geom_line() +
            labs(x = "", y = expression("Area Occupied ("~nmi^2~")")) +
            theme_bw() +
            theme(legend.title = element_blank()) 

ggsave("./figures/tanner_area_occupied.png", d95_plot,
       height = 6, width = 10)

# Write output for D95 indicator     
d95 %>%
  select(-CPUE) %>%
  pivot_wider(names_from = "CATEGORY", values_from = "d95") %>%
  write.csv("./outputs/tanner_area_occupied.csv", row.names = FALSE)



## Plot d95 vs. abund
d95_v_abund_plot <- ggplot(data = d95 %>% filter(CATEGORY != "population"), 
                          aes(x = CPUE, y = d95, group = CATEGORY, color = CATEGORY)) +
                   geom_point() +
                   # geom_line() +
                   geom_smooth(method = 'lm') +
                   labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
                   theme_bw() +
                   theme(legend.title = element_blank()) +
                   facet_wrap(~CATEGORY, scales = "free")

ggsave("./figures/tanner_area_v_abund.png", d95_v_abund_plot,
       height = 6, width = 10)


