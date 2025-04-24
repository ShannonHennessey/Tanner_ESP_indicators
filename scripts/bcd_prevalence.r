## Purpose: To calculate the visual prevalence of Bitter Crab Disease (BCD) in 
##          Tanner crab via area-expanded catch per unit effort (abundance) for
##          small crab (< 70mm), as well as for the whole EBS population.


## Read in setup
source("./scripts/setup.R")


## Subset specimen by size and disease status
# mature, with bcd
tanner_mat_bcd <- tanner
tanner_mat_bcd$specimen <- tanner_mat_bcd$specimen %>%
                           filter(SIZE >= 70,
                                  DISEASE_CODE == 2) 
# immature, with bcd
tanner_imm_bcd <- tanner
tanner_imm_bcd$specimen <- tanner_imm_bcd$specimen %>%
                           filter(SIZE < 70,
                                  DISEASE_CODE == 2) 
# mature, all
tanner_mat <- tanner
tanner_mat$specimen <- tanner_mat$specimen %>%
                       filter(SIZE >= 70) 

# immature, all
tanner_imm <- tanner
tanner_imm$specimen <- tanner_imm$specimen %>%
                       filter(SIZE < 70)

# all bcd
tanner_bcd <- tanner
tanner_bcd$specimen <- tanner_bcd$specimen %>%
                       filter(DISEASE_CODE == 2) 


## Calculate abundance for population, "immature" and "mature" subsets individually
all <- calc_bioabund(crab_data = tanner,
                     species = "TANNER",
                     spatial_level = "region") %>%
       mutate(CATEGORY = "pop")
all_bcd <- calc_bioabund(crab_data = tanner_bcd,
                         species = "TANNER",
                         spatial_level = "region") %>%
           mutate(CATEGORY = "pop_bcd")
mat_bcd <- calc_bioabund(crab_data = tanner_mat_bcd,
                         species = "TANNER",
                         spatial_level = "region") %>%
           mutate(CATEGORY = "mat_bcd")
imm_bcd <- calc_bioabund(crab_data = tanner_imm_bcd,
                         species = "TANNER",
                         spatial_level = "region") %>%
           mutate(CATEGORY = "imm_bcd")
mat <- calc_bioabund(crab_data = tanner_mat,
                    species = "TANNER",
                    spatial_level = "region") %>%
       mutate(CATEGORY = "mat")
imm <- calc_bioabund(crab_data = tanner_imm,
                     species = "TANNER",
                     spatial_level = "region") %>%
       mutate(CATEGORY = "imm")


## Combine abundance estimates, calculate percent prevalence
abund <- rbind(all, all_bcd, mat_bcd, imm_bcd, mat, imm) %>%
         select(-c("ABUNDANCE_CV","ABUNDANCE_CI", 
                   "BIOMASS_MT", "BIOMASS_MT_CV", "BIOMASS_MT_CI", 
                   "BIOMASS_LBS", "BIOMASS_LBS_CV", "BIOMASS_LBS_CI")) %>%
         pivot_wider(names_from = CATEGORY, values_from = ABUNDANCE) %>%
         mutate(POP_PREVALENCE = ifelse(pop_bcd == 0, NA, (pop_bcd/pop)*100),
                LG_PREVALENCE = ifelse(pop_bcd == 0, NA, (mat_bcd/mat)*100),
                SM_PREVALENCE = ifelse(pop_bcd == 0, NA, (imm_bcd/imm)*100)) %>%
         filter(YEAR %in% years) %>%
         select(SPECIES, YEAR, REGION, POP_PREVALENCE, LG_PREVALENCE, SM_PREVALENCE)


## Write .csv for tanner crab indicator 
write.csv(abund, "./outputs/bcd_prevalence.csv", row.names = FALSE)



## Plot time series ------------------------------------------------------------
plot_dat <- abund %>%
            pivot_longer(4:6, names_to = "CATEGORY", values_to = "PREVALENCE") %>%
            right_join(., expand_grid(YEAR = years,
                                      SPECIES = "TANNER", 
                                      REGION = "EBS",
                                      CATEGORY = unique(.$CATEGORY)))

# Combined Plot 
pop_plot <- ggplot(plot_dat, aes(x = YEAR, y = PREVALENCE, group = as.factor(CATEGORY))) +
            geom_point(aes(colour = CATEGORY), size = 3) +
            geom_line(aes(colour = CATEGORY), size = 1) +
            scale_x_continuous(limits = c(min(years), current_year),
                               breaks = seq(1980, current_year, 5)) +
            scale_color_manual(labels = c("Small prevalence", "Large prevalence", 
                                          "Population prevalence"), 
                               values = c("#E69F00", "#56B4E9", "#009E73")) +
            labs(y = "Disease Prevalence (%)", x = "Year") +
            theme_bw() +
            theme(legend.text = element_text(size = 11),
                  axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 14), 
                  axis.text.y = element_text(size = 12),
                  legend.title = element_blank())
ggsave(paste0(fig_dir, "bcd_prev.png"), pop_plot,
       height = 6, width = 10)


# Faceted plot 
labels <- function(variable,value){
  return(names[value])
}

names <- list("SM_PREVALENCE" = "Small prevalence", 
              "LG_PREVALENCE" = "Large prevalence", 
              "POP_PREVALENCE" = "Population prevalence")
facet_plot <- ggplot(plot_dat, aes(x = YEAR, y = PREVALENCE, group = as.factor(CATEGORY))) +
              geom_point(aes(colour = CATEGORY), size = 3) +
              geom_line(aes(colour = CATEGORY), size = 1) +
              scale_x_continuous(limits = c(min(years), current_year),
                                 breaks = seq(1980, current_year, 5)) +
              scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
              labs(y = "Disease Prevalence (%)", x = "Year") +
              theme_bw() +
              theme(legend.position = "none",
                    #legend.text = element_text(size = 11),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10), 
                    axis.text.y = element_text(size = 12),
                    legend.title = element_blank(),
                    panel.grid = element_blank(),
                    strip.text.x = element_text(size = 12),
                    strip.background = element_rect(fill = NA, color = NA),
                    panel.border = element_rect(colour = "black", fill = NA)) +
              facet_wrap(~CATEGORY, labeller = labels)


# Just immature crab plot
ggplot(plot_dat %>% filter(CATEGORY == "SM_PREVALENCE"), 
       aes(x = YEAR, y = PREVALENCE)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(min(years), current_year),
                     breaks = seq(1980, current_year, 5)) +
  geom_hline(aes(yintercept = mean(PREVALENCE, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(PREVALENCE, na.rm = TRUE) - sd(PREVALENCE, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(PREVALENCE, na.rm = TRUE) + sd(PREVALENCE, na.rm = TRUE)), color = "green4") +
  labs(y = "Juvenile Disease\nPrevalence (%)", x = "Year") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        legend.title = element_blank())
  
ggsave(paste0(fig_dir, "bcd_imm_prev.png"), height = 2, width = 6)

