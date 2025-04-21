## Juvenile cohort progression: quantify whether the recruitment pulses are progressing toward mature/legal sizes
## proportion of one juvenile cohort (say 45-55mm CW abundance) in year 1 to the 
## next size class (55-65mm CW) abundance in year 2 (and so on) over time to give 
## relative cohort progression -- may help flag any juvenile bottlenecks that might pop up
##
## NOTES:
## - caveats with size specific selectivity
## - example years: 2003-2005, 2008-2010 successful, 2017-2019
##   - look at differences in habitat uses between successful vs. unsuccessful 
## - another indicator to get a spatial/habitat look at when cohorts appear vs. don't?
##   - home plate closure area - does this matter?
##
## PC1: 0-25
## PC2: 25-35
## PC3: 35-50
## PC4: 50-70
## PC5: 70-85
## PC6: >85 (likely mature)


## Read in setup
source("./scripts/setup.R")


## Look at abundance distribution of sizes across years
all_crab <- calc_bioabund(crab_data = tanner,
                          species = "TANNER",
                          years = years,
                          spatial_level = "region", 
                          sex = "male",
                          bin_1mm = TRUE)

ggplot(all_crab, aes(x = SIZE_1MM, y = ABUNDANCE)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(0, 180),
                     breaks = seq(0, 180, 10)) +
  theme_bw()


## Assign pseudocohorts
cohorts <- all_crab %>%
           mutate(COHORT = case_when(SIZE_1MM < 25 ~ 1, 
                                     SIZE_1MM >= 25 & SIZE_1MM < 35 ~ 2,
                                     SIZE_1MM >= 35 & SIZE_1MM < 50 ~ 3,
                                     SIZE_1MM >= 50 & SIZE_1MM < 70 ~ 4,
                                     SIZE_1MM >= 70 & SIZE_1MM < 85 ~ 5,
                                     SIZE_1MM >= 85 ~ 6,
                                     TRUE ~ NA)) %>%
          group_by(YEAR, COHORT) %>%
          summarise(ABUNDANCE = sum(ABUNDANCE)/1e6) %>%
          group_by(YEAR) %>%
          mutate(TOTAL_ABUND = sum(ABUNDANCE),
                 PROPORTION = ABUNDANCE/TOTAL_ABUND) %>%
          right_join(., expand.grid(YEAR = years,
                                    COHORT = c(1:6))) %>%
          arrange(YEAR)


ggplot(cohorts %>% filter(!COHORT == 6), 
       aes(x = YEAR, y = ABUNDANCE, colour = as.factor(COHORT))) +
  geom_line(size = 1) +
  theme_bw()

ggplot(cohorts %>% filter(!COHORT == 6), 
       aes(x = YEAR, y = PROPORTION, colour = as.factor(COHORT))) +
  geom_line(size = 1) +
  theme_bw()


## OK so, how to calculate the t+1 etc. proportion....
progression <- cohorts %>%
               mutate()



