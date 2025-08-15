## Purpose: To calculate (1) mean female Tanner Crab size at maturity and 
##          (2) the proportion of primiparous mature female Tanner Crab 
##          with full clutches.
##
## Author: Shannon Hennessey
##
## NOTES:
## - look at differences in E166 vs. W166?
## - mature size
##   - new hardshell females only? yes?
##   - split into E and W?
## - clutch fullness -- how do we make this relate to fecundity more??
##   - lg half full female might still be more fecund than small full female
##   - fit a model that accounts for size?
## - sperm limitation? i.e. why are we seeing low clutch fullness?
##
## - maybe do all females? (hardshell...bc maybe soft haven't extruded yet)
##   - could be better signal of state of population - we recognize that there are 
##     differences in multiparous vs. primiparous females, but this is just how output
##     is in a given year...


## Read in setup
source("./scripts/setup.R")


## Calculate mean mature female size (new hardshell only)
## weighted by abundance at size
mean_size <- tanner$specimen %>%
             filter(SEX == 2,
                    CLUTCH_SIZE > 0,
                    SHELL_CONDITION == 2) %>%
             group_by(YEAR, SIZE, STATION_ID, AREA_SWEPT) %>%
             mutate(CPUE = sum(SAMPLING_FACTOR)/AREA_SWEPT) %>%
             ungroup() %>%
             group_by(YEAR) %>%
             summarize(MEAN_SIZE = weighted.mean(SIZE, w = CPUE)) %>%
             right_join(., expand.grid(YEAR = years)) %>%
             arrange(YEAR)
             


# Plot
ggplot(mean_size, aes(x = YEAR, y = MEAN_SIZE)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Mean Mature Female\nTanner Crab Size (mm)") +
  geom_hline(aes(yintercept = mean(MEAN_SIZE, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(MEAN_SIZE, na.rm = TRUE) - sd(MEAN_SIZE, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(MEAN_SIZE, na.rm = TRUE) + sd(MEAN_SIZE, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "female_SAM.png"), height = 2, width = 6)


# Save output
mean_size %>%
  rename(year = YEAR,
         female_sam = MEAN_SIZE) %>%
  write_csv("./outputs/female_SAM.csv")



## Calculate proportion of primiparous mature female Tanner Crab with full clutches ----
#**Shell 2-5, all empties!!* 
female <- tanner
female$specimen <- female$specimen %>%
                   filter(SEX == 2,
                          CLUTCH_SIZE > 0,
                          SHELL_CONDITION %in% c(2:5))

female_abund <- calc_bioabund(crab_data = female,
                                    species = "TANNER",
                                    region = "EBS",
                                    year = years,
                                    spatial_level = "region") %>%
              mutate(ALL = ABUNDANCE) %>%
              select(YEAR, ALL)


female_empty <- tanner
female_empty$specimen <- female_empty$specimen %>%
                         filter(SEX == 2,
                                CLUTCH_SIZE > 0,
                                SHELL_CONDITION %in% c(2:5),
                                EGG_CONDITION %in% c(0,3,4))

female_empty_abund <- calc_bioabund(crab_data = female_empty,
                                    species = "TANNER",
                                    region = "EBS",
                                    year = years,
                                    spatial_level = "region") %>%
                      mutate(EMPTY = ABUNDANCE) %>%
                      select(YEAR, EMPTY)

# calculate proportion
clutch_dat <- full_join(female_abund, female_empty_abund) %>%
              mutate(PROP_EMPTY = EMPTY/ALL) %>%
              right_join(., expand.grid(YEAR = years)) %>%
              arrange(YEAR)


# sperm reserves - smaller clutches when relying on reserves
# MAKE SURE TO COMMUNICATE WHY USING FULL vs. EMPTY!
# - empty doesn't seem to be hugely prevalent...does it coincide with population shifts?
# - pattern in full vs. intermediate clutches....do we see things in pop'n with that 1992 shift?

# Plot
ggplot(clutch_dat, aes(x = YEAR, y = PROP_EMPTY)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Proportion Females \nwith Empty Clutch") +
  geom_hline(aes(yintercept = mean(PROP_EMPTY, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(PROP_EMPTY, na.rm = TRUE) - sd(PROP_EMPTY, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(PROP_EMPTY, na.rm = TRUE) + sd(PROP_EMPTY, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "clutch_empty.png"), height = 2, width = 6)


# Save output
clutch_dat %>%
  select(YEAR, PROP_EMPTY) %>%
  rename(year = YEAR,
         clutch_fullness = PROP_EMPTY) %>%
  write_csv("./outputs/clutch_empty.csv")


# # calculate abundance of all primiparous mature females
# primip <- tanner
# primip$specimen <- primip$specimen %>%
#                    filter(SEX == 2,
#                           CLUTCH_SIZE > 0,
#                           SHELL_CONDITION %in% c(2))
# 
# primip_abund <- calc_bioabund(crab_data = primip,
#                               species = "TANNER",
#                               region = "EBS",
#                               year = years,
#                               spatial_level = "region") %>%
#                 mutate(PRIMIP = ABUNDANCE) %>%
#                 select(YEAR, PRIMIP)
# 
# 
# # calculate abundance of just primiparous mature females with full clutches
# primip_full <- tanner
# primip_full$specimen <- primip_full$specimen %>%
#                         filter(SEX == 2,
#                                # CLUTCH_SIZE > 0,
#                                SHELL_CONDITION %in% c(2),
#                                CLUTCH_SIZE %in% c(5:6)) #2:4 for intermediates
# 
# primip_full_abund <- calc_bioabund(crab_data = primip_full,
#                                    species = "TANNER",
#                                    region = "EBS",
#                                    year = years,
#                                    spatial_level = "region") %>%
#                      mutate(FULL_PRIMIP = ABUNDANCE) %>%
#                      select(YEAR, FULL_PRIMIP)
# 
# 
# # calculate abundance of just primiparous mature females with empty clutches
# primip_empty <- tanner
# primip_empty$specimen <- primip_empty$specimen %>%
#                          mutate() %>% # replace dead eggs with barren...think on definition
#                          filter(SEX == 2,
#                                 # CLUTCH_SIZE > 0,
#                                 SHELL_CONDITION %in% c(2),
#                                 CLUTCH_SIZE %in% c(1),
#                                 EGG_CONDITION %in% c(0,3,4))
# 
# primip_empty_abund <- calc_bioabund(crab_data = primip_empty,
#                                     species = "TANNER",
#                                     region = "EBS",
#                                     year = years,
#                                     spatial_level = "region") %>%
#                       mutate(EMPTY_PRIMIP = ABUNDANCE) %>%
#                       select(YEAR, EMPTY_PRIMIP)
# 
# # calculate proportion
# clutch_dat <- full_join(primip_full_abund, primip_abund) %>%
#               full_join(., primip_empty_abund) %>%
#               mutate(PROP_FULL = FULL_PRIMIP/PRIMIP,
#                      PROP_EMPTY = EMPTY_PRIMIP/PRIMIP) %>%
#               right_join(., expand.grid(YEAR = years)) %>%
#               arrange(YEAR) %>%
#               filter(YEAR >= 1998)
#                
# 
# 
# # Plot
# ggplot(clutch_dat, aes(x = YEAR, y = PROP_FULL)) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Year", y = "Proportion Primiparous\nFemales with Full Clutch") +
#   geom_hline(aes(yintercept = mean(PROP_FULL, na.rm = TRUE)), linetype = 5) +
#   geom_hline(aes(yintercept = mean(PROP_FULL, na.rm = TRUE) - sd(PROP_FULL, na.rm = TRUE)), color = "green4") +
#   geom_hline(aes(yintercept = mean(PROP_FULL, na.rm = TRUE) + sd(PROP_FULL, na.rm = TRUE)), color = "green4") +
#   xlim(min(years), max(years)) +
#   theme_bw()
# ggsave(paste0(fig_dir, "clutch_fullness.png"), height = 2, width = 6)
# 
# # sperm reserves - smaller clutches when relying on reserves
# # MAKE SURE TO COMMUNICATE WHY USING FULL vs. EMPTY!
# # - empty doesn't seem to be hugely prevalent...does it coincide with population shifts?
# # - pattern in full vs. intermediate clutches....do we see things in pop'n with that 1992 shift?
# 
# # Plot
# ggplot(clutch_dat, aes(x = YEAR, y = PROP_EMPTY)) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Year", y = "Proportion Primiparous\nFemales with Empty Clutch") +
#   geom_hline(aes(yintercept = mean(PROP_EMPTY, na.rm = TRUE)), linetype = 5) +
#   geom_hline(aes(yintercept = mean(PROP_EMPTY, na.rm = TRUE) - sd(PROP_EMPTY, na.rm = TRUE)), color = "green4") +
#   geom_hline(aes(yintercept = mean(PROP_EMPTY, na.rm = TRUE) + sd(PROP_EMPTY, na.rm = TRUE)), color = "green4") +
#   xlim(min(years), max(years)) +
#   theme_bw()
# ggsave(paste0(fig_dir, "clutch_empty.png"), height = 2, width = 6)
# 
# 
# # Save output
# clutch_dat %>%
#   select(YEAR, PROP_FULL) %>%
#   rename(year = YEAR,
#          clutch_fullness = PROP_FULL) %>%
#   write_csv("./outputs/clutch_fullness.csv")

