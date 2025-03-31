## Purpose: To calculate (1) Tanner Crab center of abundance (lat & lon) in EBS 
##          by size/sex category, (2) Area Occupied (D95) - area of stations 
##          that make up 95% of the cumulative Tanner CPUE, and (3) the fraction 
##          of Tanner Crab stock (abundance) in Bristol Bay 
##
## Author: Shannon Hennessey; adapted from Erin Fedewa's snow crab indicator
##
## NOTES:
## - assign male maturity by actual proportion in a given size bin, not just a cutline?
##   ...although harder to implement for years with no cutline...
##   ...but could do a global cutline for missing years? (ie fit model to all crab/all years?)
## - need to figure out what to do when missing stations in a year (see notes below for D95)
## - fraction BB: 
##    - total EBS, or just E166??
##    - not sure how to do stratification...can I just do station-level CPUE fraction?
##    - BB management district vs. BB proper?


## Read in setup
source("./scripts/setup.R")


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

## Assign maturity to specimen data; calculate CPUE
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
        summarise(COUNT = round(sum(SAMPLING_FACTOR))) %>%
        pivot_wider(names_from = CATEGORY, values_from = COUNT) %>%
        mutate(population = sum(immature_male, mature_male, immature_female, mature_female, na.rm = T)) %>%
        pivot_longer(c(6:10), names_to = "CATEGORY", values_to = "COUNT") %>%
        filter(CATEGORY != "NA") %>%
        mutate(COUNT = replace_na(COUNT, 0),
               CPUE = COUNT / AREA_SWEPT) %>%
        ungroup() 



## Compute Tanner crab center of abundance by size/sex -------------------------
# Look at # of standard hauls by year
n_haul <- tanner$haul %>% 
          group_by(YEAR) %>% 
          filter(!STATION_ID %in% corners) %>% 
          summarise(N_STATION = n())


# Calculate center of abundance
COD <- cpue %>%
       filter(!STATION_ID %in% corners) %>% # exclude corner stations
       group_by(YEAR, CATEGORY) %>%
       summarise(LAT_COD = weighted.mean(LATITUDE, w = CPUE),
                 LON_COD = weighted.mean(LONGITUDE, w = CPUE)) %>%
       right_join(., expand_grid(YEAR = years,
                                 CATEGORY = unique(cpue$CATEGORY))) %>%
       arrange(YEAR)


# Plot latitude centroid
lat_cod <- ggplot(data = COD %>% filter(CATEGORY != "population"),
                  aes(x = YEAR, y = LAT_COD, group = CATEGORY, color = CATEGORY)) +
           geom_point() +
           geom_line() +
           labs(x = "Year", y = expression(paste("Center of Abundance (", degree, "Latitude)"))) +           
           theme_bw() +
           theme(legend.title = element_blank()) 

ggsave("./figures/tanner_centroid_latitude.png", lat_cod,
       height = 6, width = 10)


# Plot longitude centroid
lon_cod <- ggplot(data = COD %>% filter(CATEGORY != "population"),
                  aes(x = YEAR, y = LON_COD, group = CATEGORY, color = CATEGORY)) +
           geom_point() +
           geom_line() +
           labs(x = "Year", y = expression(paste("Center of Abundance (", degree, "Longitude)"))) +           
           theme_bw() +
           theme(legend.title = element_blank()) 

ggsave("./figures/tanner_centroid_longitude.png", lon_cod,
       height = 6, width = 10)


# Write output for COD indicator     
COD %>%
  pivot_wider(names_from = "CATEGORY", values_from = c("LAT_COD", "LON_COD")) %>%
  write.csv("./outputs/tanner_centroid_abundance.csv", row.names = FALSE)




## Compute D95 by each size and sex category -----------------------------------
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
       filter(!STATION_ID %in% corners) %>% # exclude corner stations
       nest(data = c(-YEAR, -CATEGORY)) %>%
       mutate(d95 = purrr::map_dbl(data, f_d95_est)) %>% #apply d95 function to each element 
       unnest(cols = c(data)) %>%
       group_by(YEAR, CATEGORY) %>%
       summarise(CPUE = sum(COUNT) / sum(AREA_SWEPT), # add a column for total cpue of each group in each year
                 d95 = mean(d95)) %>% # take 'mean' just to get one value (they are all the same)
       right_join(., expand_grid(YEAR = years,
                                 CATEGORY = unique(cpue$CATEGORY))) %>%
       arrange(YEAR)

# Plot
d95_plot <- ggplot(data = d95 %>% filter(CATEGORY != "population"), 
                   aes(x = YEAR, y = d95, group = CATEGORY, color = CATEGORY)) +
            geom_point() +
            geom_line() +
            labs(x = "Year", y = expression("Area Occupied ("~nmi^2~")")) +
            theme_bw() +
            theme(legend.title = element_blank()) 

ggsave("./figures/tanner_area_occupied.png", d95_plot,
       height = 6, width = 10)


# Write output for D95 indicator     
d95 %>%
  select(-CPUE) %>%
  pivot_wider(names_from = "CATEGORY", values_from = "d95") %>%
  write.csv("./outputs/tanner_area_occupied.csv", row.names = FALSE)



test <- full_join(d95 %>% rename(year = YEAR), indicators)

# Plot d95 vs. abund
d95_v_abund_plot <- ggplot(data = d95 %>% filter(CATEGORY != "population"), 
                          aes(x = CPUE, y = d95, group = CATEGORY, color = CATEGORY)) +
                   geom_point() +
                   # geom_line() +
                   geom_smooth(method = 'lm') +
                   labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
                   theme_bw() +
                   theme(legend.title = element_blank()) +
                   facet_wrap(~CATEGORY, scales = "free")


ggplot(test %>% filter(CATEGORY != "population"), 
       aes(x = CPUE, y = d95, group = CATEGORY, color = summer_bt, label = year)) +
  geom_point() +
  # geom_line() +
  geom_smooth(method = 'lm') +
  labs(x = "CPUE", y = expression("Area Occupied ("~nmi^2~")")) +
  # lims(x = c(0,7000)) +
  geom_text(hjust = 0.5, vjust = -0.5, size = 3) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  facet_wrap(~CATEGORY, scales = "free")

ggsave("./figures/tanner_area_v_abund.png", d95_v_abund_plot,
       height = 6, width = 10)



## Calculate the fraction of Tanner stock in Bristol Bay -----------------------
# set Bristol Bay stations (136)
BB <- c("A-02","A-03","A-04","A-05","A-06","B-01","B-02","B-03","B-04","B-05","B-06",
        "B-07","B-08","C-01","C-02","C-03","C-04","C-05","C-06","C-07","C-08","C-09","D-01",
        "D-02","D-03","D-04","D-05","D-06","D-07","D-08","D-09","D-10","E-01","E-02","E-03",
        "E-04","E-05","E-06","E-07","E-08","E-09","E-10","E-11","E-12","F-01","F-02","F-03",
        "F-04","F-05","F-06","F-07","F-08","F-09","F-10","F-11","F-12","F-13","F-14","G-01",
        "G-02","G-03","G-04","G-05","G-06","G-07","G-08","G-09","G-10","G-11","G-12","G-13",
        "G-14","G-15","H-01","H-02","H-03","H-04","H-05","H-06","H-07","H-08","H-09","H-10",
        "H-11","H-12","H-13","H-14","H-15","H-16","I-01","I-02","I-03","I-04","I-05","I-06",
        "I-07","I-08","I-09","I-10","I-11","I-12","I-13","I-14","I-15","I-16","J-01","J-02",
        "J-03","J-04","J-05","J-06","J-07","J-08","J-09","J-10","J-11","J-12","J-13","J-14",
        "J-15","J-16","K-01","K-02","K-03","K-04","K-05","K-06","K-07","K-08","K-09","K-10",
        "K-11","K-12","K-13","K-14","Z-05")

## Calculate station-level CPUE for legal males
cpue_bb <- calc_cpue(crab_data = tanner,
                     species = "TANNER",
                     region = "EBS", 
                     years = years,
                     crab_category = "legal_male") %>%
          # assign stations to Bristol Bay
          mutate(STRATUM2 = ifelse(STATION_ID %in% BB, "BB", "NOT_BB")) %>%
          # remove corner stations
          filter(!STATION_ID %in% corners) %>%
          group_by(YEAR, STRATUM2) %>%
          summarise(CPUE = sum(CPUE)) %>%
          pivot_wider(names_from = STRATUM2, values_from = CPUE) %>%
          mutate(fraction_bb = BB/(BB + NOT_BB)) %>%
          right_join(., expand.grid(YEAR = years)) %>%
          arrange(YEAR)


# Plot
cpue_bb %>%
  ggplot(aes(x = YEAR, y = fraction_bb)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Fraction of Tanner Crab Stock in Bristol Bay") +
  geom_hline(aes(yintercept = mean(fraction_bb, na.rm = TRUE)), linetype = 5) +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave("./figures/fraction_bb.png")


## Save output
cpue_bb %>%
  select(YEAR, fraction_bb) %>%
  rename(year = YEAR) %>%
  write_csv("./outputs/fraction_bb.csv")



