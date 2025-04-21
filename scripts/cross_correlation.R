## Preliminary assessment of relationships...
## - need to think about lags...
## - tanner abundance as response
## - in general, lag 3 years for juv indicators, 7 for larval? could also be more though...

## Read in setup
source("./scripts/setup.R")

indicators <- read.csv("./data/BAS_indicators.csv")

recruits <- calc_bioabund(crab_data = tanner,
                          species = "TANNER",
                          years = years,
                          size_min = 70,
                          size_max = 85,
                          sex = "male",
                          spatial_level = "region") %>%
            mutate(ABUND_MIL = ABUNDANCE/1e6) %>%
            right_join(., expand.grid(YEAR = years)) %>%  
            arrange(YEAR) %>%
            mutate(LAG3 = lag(ABUND_MIL, 3),
                   LAG7 = lag(ABUND_MIL, 7))


ggplot(recruits, aes(x = YEAR, y = ABUND_MIL)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Juvenile male Tanner crab abundance (millions)") +
  geom_hline(aes(yintercept = mean(ABUND_MIL, na.rm = TRUE)), linetype = 5) +
  xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "recruits.png"), height = 4, width = 6)


# ggplot(recruits, aes(x = YEAR, y = LAG3)) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Year", y = "Juvenile male Tanner crab abundance (millions)", title = "3-year lag") +
#   geom_hline(aes(yintercept = mean(ABUND_MIL, na.rm = TRUE)), linetype = 5) +
#   xlim(min(years), max(years)) +
#   theme_bw()
# ggsave(paste0(fig_dir, "recruits_lag3.png"), height = 4, width = 6)
# 
# 
# ggplot(recruits, aes(x = YEAR, y = LAG7)) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Year", y = "Juvenile male Tanner crab abundance (millions)", title = "7-year lag") +
#   geom_hline(aes(yintercept = mean(ABUND_MIL, na.rm = TRUE)), linetype = 5) +
#   xlim(min(years), max(years)) +
#   theme_bw()
# ggsave(paste0(fig_dir, "recruits_lag7.png"), height = 4, width = 6)



## Larval indicators 
## Response: recruit abundance
## Lag: 3-5 years
# Chlorophyll-a concentration - mean_chla
# -- Wind? northeasterly favors larval advection...
# Mean summer surface temperature - summer_st
# Aleutian Low - mean_AL

## Juvenile indicators
## Response: recruit abundance
## Lag: 1
# Mean summer bottom temperature - summer_bt
# Summer cold pool extent - cp_extent
# Summer juvenile Tanner temperature occupancy - temp_occ
# DFA - temperature - dfa_temp
# Predator density - total_pred
# Daily summer pacific cod consumption of Tanner crab - pcod_consumption
# -- DFA - predation?
# Summer juvenile Tanner disease prevalence - bcs_prevalence
# -- Juvenile cohort progression? or recruitment propagation index?

## Adult indicators
## Response: recruit abundance
## Lag: 7-8 years?
# Summer benthic invertebrate density - total_invert
# Annual male Tanner size at terminal molt - male_sam
# Summer mature male Tanner area occupied (D95) - matmale_d95
# Summer mature male tanner centroid of abundance - matmale_cod_lat
# Fraction of stock in BB - fraction_bb
# -- DFA - male spatial? (^^above 3)
# Female size at maturity - female_sam
# Female clutch fullness - clutch_fullness
# -- Snow/Tanner spatial overlap?




# cross-correlation -----------------------------------------------------------
lags <- tibble(type = c(rep("larval", 3), "juvenile", rep("adult", 3)),
               lag = c(3,4,5, 1, 6,7,8))
variables <- tibble(driver = c("mean_chla","summer_st","mean_AL",
                               "summer_bt","cp_extent","temp_occ",#"dfa_temp",
                               "total_pred","pcod_consumption","bcd_prevalence",  
                               "total_invert","female_sam","clutch_fullness","male_sam",
                               "matmale_d95","matmale_cod_lat","matmale_cod_lon","fraction_bb"),
                    type = c(rep("larval", 3), rep("juvenile", 6), rep("adult", 8)), 
                    response = c(rep("recruitment", 9), NA, rep("recruitment", 2), rep(NA, 5))) #%>%
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
            rename(year = YEAR,
                   response = ABUNDANCE)

correlations <- c()

for(i in 1:nrow(variables)){

  # set up model-specific response/drivers w/ no lag
  # myresponse <- dat[, response[i]]
  driver <- indicators %>%
            select(year,
                   variables$driver[i]) %>%
            rename(driver = 2) 
  corr_df <- full_join(response, driver) %>%
             arrange(year) %>%
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

write.csv(correlations, paste0(data_dir, "cross_correlations.csv"), row.names = FALSE)

