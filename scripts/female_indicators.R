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


## Read in setup
source("./scripts/setup.R")


## Calculate mean mature female size (new hardshell only)
female_size <- calc_bioabund(crab_data = tanner,
                             species = "TANNER",
                             region = "EBS",
                             year = years,                             
                             crab_category = "mature_female",
                             shell_condition = "new_hardshell",
                             bin_1mm = TRUE)

mean_size <- female_size %>%
             group_by(YEAR) %>%
             summarize(MEAN_SIZE = weighted.mean(SIZE_1MM, w = ABUNDANCE)) %>%
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



## Calculate proportion of primiparous mature female Tanner Crab with full clutches
# calculate abundance of all primiparous mature females
primip <- tanner
primip$specimen <- primip$specimen %>%
                   filter(SEX == 2,
                          CLUTCH_SIZE > 0,
                          SHELL_CONDITION %in% c(0:2))

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
                               CLUTCH_SIZE > 0,
                               SHELL_CONDITION %in% c(0:2),
                               CLUTCH_SIZE %in% c(5:6))

primip_full_abund <- calc_bioabund(crab_data = primip_full,
                                   species = "TANNER",
                                   region = "EBS",
                                   year = years,
                                   spatial_level = "region") %>%
                     mutate(FULL_PRIMIP = ABUNDANCE) %>%
                     select(YEAR, FULL_PRIMIP)

# calculate proportion
clutch_dat <- full_join(primip_full_abund, primip_abund) %>%
              mutate(PROP_FULL = FULL_PRIMIP/PRIMIP) %>%
              right_join(., expand.grid(YEAR = years)) %>%
              arrange(YEAR)
               


# Plot
ggplot(clutch_dat, aes(x = YEAR, y = PROP_FULL)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Proportion Primiparous\nFemales with Full Clutch") +
  geom_hline(aes(yintercept = mean(PROP_FULL, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(PROP_FULL, na.rm = TRUE) - sd(PROP_FULL, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(PROP_FULL, na.rm = TRUE) + sd(PROP_FULL, na.rm = TRUE)), color = "green4") +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "clutch_fullness.png"), height = 2, width = 6)


# Save output
clutch_dat %>%
  select(YEAR, PROP_FULL) %>%
  rename(year = YEAR,
         clutch_fullness = PROP_FULL) %>%
  write_csv("./outputs/clutch_fullness.csv")

