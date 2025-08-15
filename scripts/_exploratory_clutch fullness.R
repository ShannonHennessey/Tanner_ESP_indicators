

## PRIMIPAROUS FEMALES ---------------------------------------------------------
# calculate abundance of all primiparous mature females
primip <- tanner
primip$specimen <- primip$specimen %>%
                    filter(SEX == 2,
                           CLUTCH_SIZE > 0,
                           SHELL_CONDITION %in% c(2))

primip_abund <- calc_bioabund(crab_data = primip,
                              species = "TANNER",
                              region = "EBS",
                              year = years,
                              spatial_level = "region") %>%
                mutate(PRIMIP = ABUNDANCE) %>%
                select(YEAR, PRIMIP)


# calculate abundance of just primiparous mature females with full clutches
primip_full <- tanner
primip_full$specimen <- primip_full$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(2),
         CLUTCH_SIZE %in% c(5:6))

primip_full_abund <- calc_bioabund(crab_data = primip_full,
                                   species = "TANNER",
                                   region = "EBS",
                                   year = years,
                                   spatial_level = "region") %>%
  mutate(FULL_PRIMIP = ABUNDANCE) %>%
  select(YEAR, FULL_PRIMIP)

# calculate abundance of just primiparous mature females with medium clutches
primip_med <- tanner
primip_med$specimen <- primip_med$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(2),
         CLUTCH_SIZE %in% c(2:4)) 

primip_med_abund <- calc_bioabund(crab_data = primip_med,
                                   species = "TANNER",
                                   region = "EBS",
                                   year = years,
                                   spatial_level = "region") %>%
  mutate(MED_PRIMIP = ABUNDANCE) %>%
  select(YEAR, MED_PRIMIP)


# calculate abundance of just primiparous mature females with empty clutches
primip_empty <- tanner
primip_empty$specimen <- primip_empty$specimen %>%
  mutate() %>% # replace dead eggs with barren...think on definition
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(2),
         CLUTCH_SIZE %in% c(1),
         EGG_CONDITION %in% c(0,4))

primip_empty_abund <- calc_bioabund(crab_data = primip_empty,
                                    species = "TANNER",
                                    region = "EBS",
                                    year = years,
                                    spatial_level = "region") %>%
  mutate(EMPTY_PRIMIP = ABUNDANCE) %>%
  select(YEAR, EMPTY_PRIMIP)

# calculate proportion
primip_dat <- full_join(primip_full_abund, primip_abund) %>%
  full_join(., primip_med_abund) %>%
  full_join(., primip_empty_abund) %>%
  mutate(PROP_FULL_PRIMIP = FULL_PRIMIP/PRIMIP,
         PROP_MED_PRIMIP = MED_PRIMIP/PRIMIP,
         PROP_EMPTY_PRIMIP = EMPTY_PRIMIP/PRIMIP) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  select(YEAR, PROP_FULL_PRIMIP, PROP_MED_PRIMIP, PROP_EMPTY_PRIMIP) %>%
  pivot_longer(c(2:4), names_to = "SERIES", values_to = "PROPORTION")


ggplot(primip_dat, aes(x = YEAR, y = PROPORTION, colour = SERIES)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Proportion Primiparous\nFemales") +
  xlim(min(years), max(years)) +
  theme_bw()
# ggsave(paste0(fig_dir, "clutch_fullness.png"), height = 2, width = 6)



## MULTIPAROUS FEMALES ---------------------------------------------------------
# calculate abundance of all primiparous mature females
primip <- tanner
primip$specimen <- primip$specimen %>%
  filter(SEX == 2,
         CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(3:5))

primip_abund <- calc_bioabund(crab_data = primip,
                              species = "TANNER",
                              region = "EBS",
                              year = years,
                              spatial_level = "region") %>%
  mutate(PRIMIP = ABUNDANCE) %>%
  select(YEAR, PRIMIP)


# calculate abundance of just primiparous mature females with full clutches
primip_full <- tanner
primip_full$specimen <- primip_full$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(3:5),
         CLUTCH_SIZE %in% c(5:6))

primip_full_abund <- calc_bioabund(crab_data = primip_full,
                                   species = "TANNER",
                                   region = "EBS",
                                   year = years,
                                   spatial_level = "region") %>%
  mutate(FULL_PRIMIP = ABUNDANCE) %>%
  select(YEAR, FULL_PRIMIP)

# calculate abundance of just primiparous mature females with medium clutches
primip_med <- tanner
primip_med$specimen <- primip_med$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(3:5),
         CLUTCH_SIZE %in% c(2:4)) 

primip_med_abund <- calc_bioabund(crab_data = primip_med,
                                  species = "TANNER",
                                  region = "EBS",
                                  year = years,
                                  spatial_level = "region") %>%
  mutate(MED_PRIMIP = ABUNDANCE) %>%
  select(YEAR, MED_PRIMIP)


# calculate abundance of just primiparous mature females with empty clutches
primip_empty <- tanner
primip_empty$specimen <- primip_empty$specimen %>%
  mutate() %>% # replace dead eggs with barren...think on definition
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(3:5),
         CLUTCH_SIZE %in% c(1),
         EGG_CONDITION %in% c(0,4))

primip_empty_abund <- calc_bioabund(crab_data = primip_empty,
                                    species = "TANNER",
                                    region = "EBS",
                                    year = years,
                                    spatial_level = "region") %>%
  mutate(EMPTY_PRIMIP = ABUNDANCE) %>%
  select(YEAR, EMPTY_PRIMIP)

# calculate proportion
multip_dat <- full_join(primip_full_abund, primip_abund) %>%
  full_join(., primip_med_abund) %>%
  full_join(., primip_empty_abund) %>%
  mutate(PROP_FULL_MULTIP = FULL_PRIMIP/PRIMIP,
         PROP_MED_MULTIP = MED_PRIMIP/PRIMIP,
         PROP_EMPTY_MULTIP = EMPTY_PRIMIP/PRIMIP) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  select(YEAR, PROP_FULL_MULTIP, PROP_MED_MULTIP, PROP_EMPTY_MULTIP) %>%
  pivot_longer(c(2:4), names_to = "SERIES", values_to = "PROPORTION")


ggplot(multip_dat, aes(x = YEAR, y = PROPORTION, colour = SERIES)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Proportion Multiparous\nFemales") +
  xlim(min(years), max(years)) +
  theme_bw()
# ggsave(paste0(fig_dir, "clutch_fullness.png"), height = 2, width = 6)


## ALL FEMALES ---------------------------------------------------------
# calculate abundance of all primiparous mature females
primip <- tanner
primip$specimen <- primip$specimen %>%
  filter(SEX == 2,
         CLUTCH_SIZE > 0)

primip_abund <- calc_bioabund(crab_data = primip,
                              species = "TANNER",
                              region = "EBS",
                              year = years,
                              spatial_level = "region") %>%
  mutate(PRIMIP = ABUNDANCE) %>%
  select(YEAR, PRIMIP)


# calculate abundance of just primiparous mature females with full clutches
primip_full <- tanner
primip_full$specimen <- primip_full$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         # SHELL_CONDITION %in% c(2),
         CLUTCH_SIZE %in% c(5:6))

primip_full_abund <- calc_bioabund(crab_data = primip_full,
                                   species = "TANNER",
                                   region = "EBS",
                                   year = years,
                                   spatial_level = "region") %>%
  mutate(FULL_PRIMIP = ABUNDANCE) %>%
  select(YEAR, FULL_PRIMIP)

# calculate abundance of just primiparous mature females with medium clutches
primip_med <- tanner
primip_med$specimen <- primip_med$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         # SHELL_CONDITION %in% c(2),
         CLUTCH_SIZE %in% c(2:4)) 

primip_med_abund <- calc_bioabund(crab_data = primip_med,
                                  species = "TANNER",
                                  region = "EBS",
                                  year = years,
                                  spatial_level = "region") %>%
  mutate(MED_PRIMIP = ABUNDANCE) %>%
  select(YEAR, MED_PRIMIP)


# calculate abundance of just primiparous mature females with empty clutches
primip_empty <- tanner
primip_empty$specimen <- primip_empty$specimen %>%
  mutate() %>% # replace dead eggs with barren...think on definition
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         # SHELL_CONDITION %in% c(2),
         CLUTCH_SIZE %in% c(1),
         EGG_CONDITION %in% c(0,4))

primip_empty_abund <- calc_bioabund(crab_data = primip_empty,
                                    species = "TANNER",
                                    region = "EBS",
                                    year = years,
                                    spatial_level = "region") %>%
  mutate(EMPTY_PRIMIP = ABUNDANCE) %>%
  select(YEAR, EMPTY_PRIMIP)

# calculate proportion
female_dat <- full_join(primip_full_abund, primip_abund) %>%
  full_join(., primip_med_abund) %>%
  full_join(., primip_empty_abund) %>%
  mutate(PROP_FULL_ALL = FULL_PRIMIP/PRIMIP,
         PROP_MED_ALL = MED_PRIMIP/PRIMIP,
         PROP_EMPTY_ALL = EMPTY_PRIMIP/PRIMIP) %>%
  right_join(., expand.grid(YEAR = years)) %>%
  arrange(YEAR) %>%
  select(YEAR, PROP_FULL_ALL, PROP_MED_ALL, PROP_EMPTY_ALL) %>%
  pivot_longer(c(2:4), names_to = "SERIES", values_to = "PROPORTION")


ggplot(female_dat, aes(x = YEAR, y = PROPORTION, colour = SERIES)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Proportion Females") +
  xlim(min(years), max(years)) +
  theme_bw()
# ggsave(paste0(fig_dir, "clutch_fullness.png"), height = 2, width = 6)

## EXPLORATION -----------------------------------------------------------------

dat <- rbind(primip_dat, multip_dat, female_dat) %>%
       pivot_wider(names_from = SERIES, values_from = PROPORTION)


# mean clutch fullness
# primip
primip <- tanner
primip <- primip$specimen %>%
          filter(SEX == 2,
                 CLUTCH_SIZE > 0,
                 SHELL_CONDITION %in% c(2)) %>%
          group_by(YEAR) %>%
          summarise(MEAN_PRIMIP = weighted.mean(CLUTCH_SIZE, SAMPLING_FACTOR))

# multip
multip <- tanner
multip <- multip$specimen %>%
  filter(SEX == 2,
         CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(3:5)) %>%
  group_by(YEAR) %>%
  summarise(MEAN_MULTIP = weighted.mean(CLUTCH_SIZE, SAMPLING_FACTOR))


# all
all <- tanner
all <- all$specimen %>%
  filter(SEX == 2,
         CLUTCH_SIZE > 0) %>%
  group_by(YEAR) %>%
  summarise(MEAN_ALL = weighted.mean(CLUTCH_SIZE, SAMPLING_FACTOR))


ggplot(full_join(primip, multip) %>%
         full_join(., all) %>%
         pivot_longer(c(2:4), names_to = "SERIES", values_to = "PROPORTION"), 
       aes(x = YEAR, y = PROPORTION, colour = SERIES)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Mean clutch size") +
  xlim(min(years), max(years)) +
  theme_bw()



dat <- dat %>% 
  full_join(., primip) %>%
  full_join(., multip) %>%
  full_join(., all)


# cross-correlation -----------------------------------------------------------
lags <- tibble(type = c(rep("larval", 3), "juvenile", rep("adult", 3)),
               lag = c(3,4,5, 1, 6,7,8))
variables <- tibble(driver = c("PROP_FULL_PRIMIP", "PROP_MED_PRIMIP", "PROP_EMPTY_PRIMIP",
                               "PROP_FULL_MULTIP", "PROP_MED_MULTIP", "PROP_EMPTY_MULTIP",
                               "PROP_FULL_ALL", "PROP_MED_ALL", "PROP_EMPTY_ALL",
                               "MEAN_PRIMIP", "MEAN_MULTIP", "MEAN_ALL"),
                    type = c(rep("adult", 12)), 
                    response = c(rep("recruitment", 12))) #%>%
# right_join(., lags)

response <- calc_bioabund(crab_data = tanner,
                          species = "TANNER",
                          spatial_level = "region",
                          years = years,
                          sex = "male",
                          size_min = 70, 
                          size_max = 85) %>%
            select(YEAR, ABUNDANCE) %>%
            right_join(., expand.grid(YEAR = years)) %>%
            arrange(YEAR) %>%
            rename(response = ABUNDANCE)

correlations <- c()

for(i in 1:nrow(variables)){
  
  # set up model-specific response/drivers w/ no lag
  # myresponse <- dat[, response[i]]
  driver <- dat %>%
            select(YEAR,
                   variables$driver[i]) %>%
            rename(driver = 2) 
  corr_df <- full_join(response, driver) %>%
             arrange(YEAR) %>%
             rename(year = YEAR) %>%
             filter(!(is.na(response) & is.na(driver) & year < 1982))
  
  # calculate cross-correlations
  cross_corr <- ccf(corr_df$driver, corr_df$response, 
                    na.action = na.pass, lag.max = 10)
  
  acf <- cross_corr$acf[which(abs(cross_corr$acf) == max(abs(cross_corr$acf)))]
  lag <- ifelse(length(cross_corr$lag[which(cross_corr$acf == acf)]) == 0, NA,
                cross_corr$lag[which(cross_corr$acf == acf)])
  
  # estimate significance
  ci <- qnorm((1 + 0.95)/2)/sqrt(cross_corr$n.used)
  sig <- ifelse(abs(acf) > ci, 1, 0)
  
  
  # extract relevant model info/values for comparison
  temp <- c(driver = variables$driver[i],
            response = "recruitment", #variables$response[i],
            lag = ifelse(sig == 0, NA, lag),
            acf = acf)
  
  correlations <- rbind(correlations, temp)
}

# reshape data frame
correlations <- as.data.frame(correlations)
# x(driver)t+k(lag) = y(response/recruitment)t 


# ------------

primip <- tanner
primip$specimen <- primip$specimen %>%
  filter(SEX == 2,
         # CLUTCH_SIZE > 0,
         SHELL_CONDITION %in% c(2))

primip_abund <- calc_bioabund(crab_data = primip,
                              species = "TANNER",
                              region = "EBS",
                              year = years,
                              spatial_level = "region",
                              clutch_size = "all_categories") %>%
                select(YEAR, CLUTCH_TEXT, ABUNDANCE) %>%
                filter(!CLUTCH_TEXT == "immature") %>%
                group_by(YEAR) %>%
                mutate(ALL = sum(ABUNDANCE),
                       PROPORTION = ABUNDANCE/ALL) 

ggplot(primip_abund, aes(x = YEAR, y = PROPORTION, colour = CLUTCH_TEXT)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Proportion Females") +
  xlim(min(years), max(years)) +
  theme_bw()

