## Purpose: To calculate tanner crab CPUE-weighted mean temperature of occupancy 
##
## NOTES:
## - need to impute missing station temperatures for early years 
## - do we think there's a behavioral shift with maturity?? Or is it more size stratified?
##   - could see the argument for maturity bc mating overlaps etc...but size-dependent predation...temp differences?


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
        summarise(COUNT = round(sum(SAMPLING_FACTOR))) %>%
        pivot_wider(names_from = CATEGORY, values_from = COUNT) %>%
        mutate(population = sum(immature_male, mature_male, immature_female, mature_female, na.rm = T)) %>%
        pivot_longer(c(6:10), names_to = "CATEGORY", values_to = "COUNT") %>%
        filter(CATEGORY != "NA") %>%
        mutate(COUNT = replace_na(COUNT, 0),
               CPUE = COUNT / AREA_SWEPT) %>%
        ungroup() 



## Calculate Temperature of Occupancy
# Ignoring NA's in bottom temp data -- missing values should be imputed in future iterations 
temp_occ <- cpue %>%
            # join in haul-level temperature data
            left_join(., tanner$haul %>% select(YEAR, GEAR_TEMPERATURE, STATION_ID)) %>%
            # calculate yearly mean bottom temperature
            group_by(YEAR) %>%
            mutate(MEAN_BT = mean(GEAR_TEMPERATURE, na.rm = T)) %>%
            ungroup() %>%
            # subset to just immature
            filter(CATEGORY %in% c("immature_male", "immature_female")) %>%
            # calculate mean bottom temperature weighted by CPUE
            group_by(YEAR, MEAN_BT) %>% # CATEGORY
            summarise(TEMP_OCC = weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T)) %>%
            right_join(., expand_grid(YEAR = years)) %>%
            arrange(YEAR) #%>%
            # mutate(DIFF = TEMP_OCC - MEAN_BT) # generally at warmer temps than mean bottom temp
            
  

## Plot
ggplot(data = temp_occ,
       aes(x = YEAR, y = TEMP_OCC)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mean(TEMP_OCC, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(TEMP_OCC, na.rm = TRUE) - sd(TEMP_OCC, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(TEMP_OCC, na.rm = TRUE) + sd(TEMP_OCC, na.rm = TRUE)), color = "green4") +
  labs(x = "Year", y = "Temperature\nOccupied (C)") + #expression(paste("Temperature Occupied ", "(", degree, "C)"))            
  theme_bw() +
  theme(legend.title = element_blank()) 
ggsave(paste0(fig_dir, "tanner_temp_occupied.png"), height = 2, width = 6)


## Write output for Temp Occupancy indicator     
temp_occ %>%
  select(-MEAN_BT) %>%
  rename(year = YEAR,
         temp_occ = TEMP_OCC) %>%
  # pivot_wider(names_from = "CATEGORY", values_from = "TEMP_OCC") %>%
  write.csv("./outputs/tanner_temp_occupied.csv", row.names = FALSE)


