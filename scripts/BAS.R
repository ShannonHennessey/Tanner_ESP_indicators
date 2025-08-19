## Project Name: ECOSYSTEM AND SOCIOECONOMIC PROFILES - Eastern Bering Sea Tanner Crab
##
## Creator: Dr. Curry James Cunningham, UAF, CFOS
## With additions from S. Hennessey and E. Fedewa, NOAA AFSC
## Date: September 2025
##
## Purpose: To evaluate linkages between recruitment and a standard set of atmospheric, 
##          oceanographic, and biological indicators of ecosystem status for the 
##          Ecosystem and Socioeconomic Profiles (ESPs).
##
## BAS RESPONSE VARIABLE
## Calculate abundance of immature male Tanner crab 70-85mm as response for BAS analysis
## Size range selected using BSFRF selectivity curves and Donaldson et al. 1981 size at 
## age estimates (~4-5 years post settlement for GOA Tanner, assuming EBS Tanner grow slower so 5-6 years?)
## Develop a female response variable too for additional model runs?


## Setup -----------------------------------------------------------------------
# Load packages
library(tidyverse)
library(crabpack)
library(ggridges)

require(corrplot)
require(cowplot)
require(viridis)
require(ggthemes)
require(BAS)
require(readxl)
require(gbm)

# Set data directory
data_dir <- "Y:/KOD_Research/Hennessey/Tanner_ESP/data/"
fig_dir <- "Y:/KOD_Research/Hennessey/Tanner_ESP/figures/BAS/"

# Set years
current_year <- 2025
years <- 1982:current_year


# Pull Tanner specimen data
tanner <- get_specimen_data(species = "TANNER",
                            region = "EBS",
                            channel = "KOD")

# Calculate abundance for response variable
response <- calc_bioabund(crab_data = tanner,
                          species = "TANNER",
                          spatial_level = "region",
                          years = years,
                          sex = "male",
                          size_min = 70, 
                          size_max = 85) %>%
           select(YEAR, ABUNDANCE) %>%
           mutate(ABUNDANCE = ABUNDANCE/1e6) %>%
           right_join(., expand.grid(YEAR = years)) %>%
           arrange(YEAR) %>%
           rename(year = YEAR,
                  imm_survey_abund = ABUNDANCE)

# Plot response to check
ggplot(response, aes(x = year, y = imm_survey_abund)) +
  geom_point() +
  geom_line() +
  theme_bw()

# Write csv as output (abundance in millions of crab)
write.csv(response, file = ("./outputs/BAS_response.csv"), row.names=FALSE)


# Pull indicator data
indicators <- read.csv("./data/BAS_indicators.csv") %>%
              mutate(total_invert = ifelse(year >= 1988, total_invert, NA))

# Join indicator and response data
dat_tanner <- indicators %>% 
              left_join(response)



## Run Model -------------------------------------------------------------------
## Using design-based BT survey estimate for male recruitment as response

# Look at temporal coverage of indicators 
dat_tanner %>%
  select(!imm_survey_abund) %>%
  pivot_longer(c(2:(ncol(dat_tanner)-1)), names_to = "indicator", values_to = "value") %>%
  ggplot(aes(year, indicator, size = value)) +
  geom_point() +
  theme_bw()
# Drop chla below, as they will constrain BAS due to short timeseries
# We'll also drop any spatial distribution indicators, as these are not drivers of recruitment 


# Set up lags
variables <- tibble(indicator = c("summer_st","wind_along_shelf","wind_cross_shelf","NPI","mean_AL","ice_avg",
                                  "mean_chla","temp_occ","total_pred","pcod_consumption","bcd_prevalence",  
                                  "total_invert","female_sam","clutch_fullness",
                                  "male_sam","matmale_d95","matmale_cod_lon"), 
                    type = c(rep("larval", 5), rep("juvenile", 7), rep("adult", 5)), 
                    response = c(rep("recruitment", 12), rep(NA, 5)),
                    lag = c(5,5,5,5,5,2, 3,1,3,3,3,2, 7,7, rep(NA, 3))) %>%
            filter(!indicator %in% c("female_sam","clutch_fullness",
                                    "male_sam","matmale_d95","matmale_cod_lon"))

# Set iteration identifier for saving outputs
iter <- "sept15"


## Assign lags for indicators - see metadata file in repo for rationales for lags
dat_tanner_bas <- dat_tanner %>%
                  select(-imm_survey_abund) %>%
                  pivot_longer(c(2:ncol(.)), names_to = "indicator", values_to = "value") %>%
                  filter(year >= 1982,
                         indicator %in% variables$indicator) %>%
                  left_join(variables) %>%
                  filter(!is.na(lag)) %>%
                  group_by(indicator) %>%
                  nest() %>%
                  mutate(data = purrr::map(data, function(data){
                    n_lag <- as.numeric(unique(data$lag))
                    x <- data %>%
                         mutate(lagged = lag(value, n = n_lag, order_by = year))
                    return(x)})) %>%
                  unnest(cols = c(data)) %>%
                  select(indicator, year, lagged) %>%
                  pivot_wider(names_from = "indicator", values_from = "lagged") %>%
                  left_join(response)

# Plot again and look at temporal coverage with lags incorporated
# with so many large ELH lags, we're going to lose early years in the timeseries
dat_tanner_bas %>%
  select(-imm_survey_abund) %>%
  pivot_longer(c(2:ncol(.)), names_to = "indicator", values_to = "value") %>%
  ggplot(aes(x = year, y = indicator, size = value)) +
  geom_point(na.rm = T) +
  theme_bw()


# # Plot timeseries with lagged covariates
# dat_tanner_bas %>%
#   pivot_longer(c(2:ncol(.)), names_to = "indicator", values_to = "value") %>%
#   ggplot(aes(year, value)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~indicator, scales = "free_y") +
#   theme_bw()
# 
# # Look at distributions of potentially problematic covariates
# hist(dat_tanner_bas$imm_survey_abund) # skew left
# hist(dat_tanner_bas$bcd_prevalence) # skew left
# hist(dat_tanner_bas$mean_AL)
# hist(dat_tanner_bas$NPI)
# hist(dat_tanner_bas$wind_along_shelf)
# hist(dat_tanner_bas$wind_cross_shelf)
# hist(dat_tanner_bas$mean_chla)
# hist(dat_tanner_bas$pcod_consumption) # skew left
# hist(dat_tanner_bas$summer_st)
# hist(dat_tanner_bas$ice_avg)
# hist(dat_tanner_bas$temp_occ)
# hist(dat_tanner_bas$total_invert)
# hist(dat_tanner_bas$total_pred)
# hist(dat_tanner_bas$female_sam)
# hist(dat_tanner_bas$clutch_fullness)
# 
# #Determine Covariates
# if(model == "BAS_Sept_2025") {
#   covars <- names(dat_tanner_bas %>% select(-year, -imm_survey_abund))
# }
# 
# 
# 
# # # Calculate Log Recruitment
# # dat_tanner_bas <- dat_tanner_bas %>%
# #                   mutate(ln_abund = log(imm_survey_abund),
# #                          pcod_consumption = log(pcod_consumption),
# #                          bcd_prevalence = log(bcd_prevalence))
# 
# # Log transform skewed predictors 
# hist(log(dat_tanner_bas$pcod_consumption))
# hist(log(dat_tanner_bas$bcd_prevalence))
# 
# # if(model == "BAS_Sep_2025") {
# #   dat_tanner_bas$consump_lag <- log(dat_tanner_bas$consump_lag)
# #   dat_tanner_bas$bcs_lag <- log(dat_tanner_bas$bcs_lag)
# # }
# 
# 
# dat_tanner_bas <- dat_tanner_bas %>%
#                   mutate(ln_abund = log(imm_survey_abund),
#                          pcod_consumption = log(pcod_consumption),
#                          bcd_prevalence = log(bcd_prevalence))
# 
# 
# 

## Final data shaping ----------------------------------------------------------
dat_tanner_bas <- dat_tanner %>%
                  # remove unused indicators
                  select(-imm_survey_abund, -mean_chla, -mean_AL, -NPI, -wind_cross_shelf, -ice_avg) %>% 
                  pivot_longer(c(2:ncol(.)), names_to = "indicator", values_to = "value") %>%
                  filter(year >= 2000, 
                         indicator %in% variables$indicator) %>%
                  left_join(variables) %>%
                  filter(!is.na(lag)) %>%
                  group_by(indicator) %>%
                  nest() %>%
                  mutate(data = purrr::map(data, function(data){
                    n_lag <- as.numeric(unique(data$lag))
                    x <- data %>% 
                         mutate(lagged = lag(value, n = n_lag, order_by = year))
                    return(x)})) %>%
                  unnest(cols = c(data)) %>%
                  select(indicator, year, lagged) %>%
                  pivot_wider(names_from = "indicator", values_from = "lagged") %>%
                  left_join(response) %>% 
                  # add log recruitment
                  mutate(ln_abund = log(imm_survey_abund))


# Determine covariates
covars <- names(dat_tanner_bas %>% select(-year, -imm_survey_abund, -ln_abund))



## Standardize covariates
# Plot
covar.list <- dat_tanner_bas %>% 
              select(-imm_survey_abund, -ln_abund) %>% 
              gather(key = type, value = value, -year) 


# Z-score predictors that are bounded at zero
dat.4 <- dat_tanner_bas
c <- 1
for(c in 1:length(covars)) {
  dat.4[[covars[c]]] <- (dat.4[[covars[c]]] - mean(dat.4[[covars[c]]], na.rm=TRUE)) / sd(dat.4[[covars[c]]], na.rm = TRUE)
}

# Make sure all predictors are correctly z-scored
apply(dat.4, 2, mean, na.rm = TRUE)
apply(dat.4, 2, sd, na.rm = TRUE)


# Subset Data for Fitting 
dat.fit <- dat.4 %>% select(-imm_survey_abund, ln_abund)
dat.fit.list <- dat.fit %>% gather(key = 'var', value = 'value', -year)


# Correlation Plot
covar.mtx <- dat.fit %>% 
             select(-year, -ln_abund)

corr.mtx <- cor(covar.mtx, use = "na.or.complete")
png(paste0(fig_dir, "BAS_covariate_correlation_", iter, ".png"), 
    height = 12, width = 12, units = 'in', res = 300)
corrplot::corrplot(corr.mtx, method = "number")
dev.off()

# Check number of pairs with correlation >0.6
sum((corr.mtx > 0.6 & corr.mtx < 1.0) & (corr.mtx < -0.6))



## Fit models ------------------------------------------------------------------
# Remove year, rename variables 
dat.temp <- dat.fit %>% 
            select(-year) %>%
            rename("Surface Temperature" = summer_st, 
                   # "ALBSA" = mean_AL, #Aleutian Low - Beaufort Sea Anticyclone
                   # "NPI" = NPI,
                   # "Sea Ice Extent" = ice_avg,
                   "Along-Shelf Wind" = wind_along_shelf,
                   # "Cross-Shelf Wind" = wind_cross_shelf,
                   # "Chlorophyll-a Concentration" = mean_chla,
                   "Occupied Temperature" = temp_occ,
                   "Predator Density" = total_pred, 
                   "Disease Prevalance" = bcd_prevalence,
                   "Benthic Prey" = total_invert,
                   # "Female Size at Maturity" = female_sam,
                   # "Clutch Failure" = clutch_fullness,
                   "Pacific Cod Consumption" = pcod_consumption)

# Bayesian model melection
bas.lm <-  bas.lm(ln_abund ~ ., 
                  data = dat.temp,
                  modelprior = uniform(), 
                  initprobs = "Uniform",
                  method = 'BAS', 
                  MCMC.iterations = 1e5, 
                  thin = 10)
summary(bas.lm)

# Diagnostic plots
plot(bas.lm, which = 4, ask = FALSE, caption = "", sub.caption = "")
plot(coef(bas.lm), ask = FALSE)
plot(confint(coef(bas.lm, level = 0.95)))
plot(bas.lm, which = 4)


## Plot model predictions vs. observed -----------------------------------------
pdf(paste0(fig_dir, "BAS_model_fit_", iter, ".pdf"), height = 5, width = 10)
par(oma = c(1,1,1,1), mar = c(4,4,1,1), mfrow = c(1,2))
pred.bas <- predict(bas.lm, estimator = "BMA")

# Omit NAs
dat.temp.na.omit <- na.omit(dat.fit)

plot(x = dat.temp.na.omit$ln_abund, y = pred.bas$Ybma,
     xlab = "Observed ln(Recruitment)", ylab = "Predicted ln(Recruitment)", 
     pch = 21, bg = rgb(1,0,0,alpha = 0.5), main = "")

# Title
mtext(paste("Tanner Crab", model), side = 3, outer = TRUE, font = 2)
abline(a = 0, b = 1, col = rgb(0,0,1,alpha = 0.5), lwd = 3)

# Timeseries
plot(x = dat.temp.na.omit$year, y = dat.temp.na.omit$ln_abund,
     xlab = "Year", ylab = "ln(Recruitment)", type = 'l', 
     col = rgb(1,0,0,alpha = 0.5), main = "")
grid(lty = 3, col = 'dark gray')
points(x = dat.temp.na.omit$year, y = dat.temp.na.omit$ln_abund,
       pch = 21, bg = rgb(1,0,0,alpha = 0.5))
lines(x = dat.temp.na.omit$year, y = pred.bas$Ybma, 
      lwd = 3, col = rgb(0,0,1, alpha = 0.5))
legend('topleft', legend = c("Observed", "Predicted"), 
       lty = 1, col = c(rgb(1,0,0, alpha = 0.5), rgb(0,0,1, alpha = 0.5)),
       bg = "white")

dev.off()



## Plot inclusion probabilities ------------------------------------------------
inc.probs <- summary(bas.lm)[2:ncol(dat.temp), 1]

bas.names <- coef(bas.lm)$namesx
inc.probs <- coef(bas.lm)$probne0
post.mean <- coef(bas.lm)$postmean
post.sd <- coef(bas.lm)$postsd
low.95 <- confint(coef(bas.lm))[,1]
up.95 <- confint(coef(bas.lm))[,2]


# Make final plot of covariates
par(mfrow = c(1, 2), mar = c(4,1,2,1), oma = c(0,10,1,1))
plot.df <- data.frame(bas.names, inc.probs, post.mean, post.sd, low.95, up.95)

# Parameter estimates/CI
g.b <- ggplot(filter(plot.df, bas.names != 'Intercept'),
              aes(x = bas.names, post.mean, fill = 'royalblue4')) +
       theme_bw() +
       geom_errorbar(aes(ymin = low.95, ymax = up.95), width = 0.25) +
       geom_point(pch = 21, fill = 'royalblue4', size = 3) +
       geom_hline(yintercept = 0, col = 'red', alpha = 0.5) +
       ylab('Effect') +
       xlab('Covariate') +
       coord_flip() +
       theme(legend.position = 'none')

# Inclusion probability
g2.b <- ggplot(filter(plot.df, bas.names != 'Intercept'),
               aes(x = bas.names, y = inc.probs, fill = inc.probs)) +
        theme_bw() +
        geom_bar(stat = 'identity', color = 'black') +
        ylab('Inclusion\nProbability') +
        scale_y_continuous(limits = c(0, 1)) +
        geom_hline(yintercept = c(0, 1)) +
        geom_hline(yintercept = 0.5, col = 'black', linetype = 5, alpha = 0.5) +
        theme(legend.position = 'none', axis.text.y = element_blank(), 
              axis.title.y = element_blank()) +
        coord_flip() +
        scale_fill_continuous_tableau()


# Bring figures together and save
g3.b <- plot_grid(g.b, g2.b, nrow = 1, ncol = 2, rel_widths = c(3, 1), align = 'h')
g3.b
ggsave(file = paste0(fig_dir, "BAS_", iter, ".png"), plot = g3.b, height = 5, width = 8, 
       units = 'in', dpi = 500)



