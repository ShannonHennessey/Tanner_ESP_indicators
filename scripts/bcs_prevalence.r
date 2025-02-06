## Purpose: To calculate the visual prevalence of Bitter Crab Syndrome (BCS) in 
##          Tanner crab via area-expanded catch per unit effort (abundance) for
##          mature crab and immature crab (based on morphometric maturity size
##          estimates), as well as for the whole EBS population.


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
            rename(SAM = B_EST, 
                   STD_ERR = B_SE) %>%
            right_join(., expand_grid(YEAR = c(1975:2024),
                                      SPECIES = "TANNER", 
                                      REGION = "EBS",
                                      DISTRICT = c("ALL", "E166", "W166"))) %>%
            mutate(SAM = case_when(DISTRICT == "ALL" & is.na(SAM) ~ 103.5, 
                                   DISTRICT == "E166" & is.na(SAM) ~ 110, 
                                   DISTRICT == "W166" & is.na(SAM) ~ 99, 
                                   TRUE ~ SAM))
  
## Assign maturity to specimen data
tanner$specimen <- tanner$specimen %>% 
                   left_join(., mat_size) %>%
                   mutate(CATEGORY = case_when((SEX == 1 & SIZE >= SAM) | (SEX == 2 & CLUTCH_SIZE >= 1) ~ "mature",
                                               (SEX == 1 & SIZE< SAM) | (SEX == 2 & CLUTCH_SIZE == 0) ~ "immature",
                                               TRUE ~ NA)) 

## Subset specimen by maturity and disease status
# mature, with BCS
tanner_mat_bcs <- tanner
tanner_mat_bcs$specimen <- tanner_mat_bcs$specimen %>%
                           filter(DISEASE_CODE == 2, 
                                  CATEGORY == "mature")
# immature, with BCS
tanner_imm_bcs <- tanner
tanner_imm_bcs$specimen <- tanner_imm_bcs$specimen %>%
                           filter(DISEASE_CODE == 2, 
                                  CATEGORY == "immature")
# mature, all
tanner_mat <- tanner
tanner_mat$specimen <- tanner_mat$specimen %>%
                       filter(CATEGORY == "mature") 
# immature, all
tanner_imm <- tanner
tanner_imm$specimen <- tanner_imm$specimen %>%
                       filter(CATEGORY == "immature")
# all BCS
tanner_bcs <- tanner
tanner_bcs$specimen <- tanner_bcs$specimen %>%
                       filter(DISEASE_CODE == 2) 


## Calculate abundance for population, "immature" and "mature" subsets individually
all <- calc_bioabund(crab_data = tanner,
                     species = "TANNER",
                     spatial_level = "region") %>%
       mutate(CATEGORY = "pop")
all_bcs <- calc_bioabund(crab_data = tanner_bcs,
                         species = "TANNER",
                         spatial_level = "region") %>%
           mutate(CATEGORY = "pop_bcs")
mat_bcs <- calc_bioabund(crab_data = tanner_mat_bcs,
                         species = "TANNER",
                         spatial_level = "region") %>%
           mutate(CATEGORY = "mat_bcs")
imm_bcs <- calc_bioabund(crab_data = tanner_imm_bcs,
                         species = "TANNER",
                         spatial_level = "region") %>%
           mutate(CATEGORY = "imm_bcs")
mat <- calc_bioabund(crab_data = tanner_mat,
                    species = "TANNER",
                    spatial_level = "region") %>%
       mutate(CATEGORY = "mat")
imm <- calc_bioabund(crab_data = tanner_imm,
                     species = "TANNER",
                     spatial_level = "region") %>%
       mutate(CATEGORY = "imm")


## Combine abundance estimates, calculate percent prevalence
abund2 <- rbind(all, all_bcs, mat_bcs, imm_bcs, mat, imm) %>%
         select(-c("ABUNDANCE_CV","ABUNDANCE_CI", 
                   "BIOMASS_MT", "BIOMASS_MT_CV", "BIOMASS_MT_CI", 
                   "BIOMASS_LBS", "BIOMASS_LBS_CV", "BIOMASS_LBS_CI")) %>%
         pivot_wider(names_from = CATEGORY, values_from = ABUNDANCE) %>%
         mutate(POP_PREVALENCE = (pop_bcs / pop)*100,
                MAT_PREVALENCE = (mat_bcs / mat)*100,
                IMM_PREVALENCE = (imm_bcs / imm)*100) %>%
         filter(YEAR >= 1989) %>%
         select(SPECIES, YEAR, REGION, POP_PREVALENCE, MAT_PREVALENCE, IMM_PREVALENCE)


## Write .csv for tanner crab indicator 
write.csv(abund, "./outputs/bcs_prevalence.csv", row.names = FALSE)



## Plot time series ------------------------------------------------------------
current_year <- 2024

plot_dat <- abund %>%
            pivot_longer(4:6, names_to = "CATEGORY", values_to = "PREVALENCE") %>%
            right_join(., expand_grid(YEAR = c(1989:current_year),
                                      SPECIES = "TANNER", 
                                      REGION = "EBS",
                                      CATEGORY = unique(.$CATEGORY)))

# Combined Plot 
pop_plot <- ggplot(plot_dat, aes(x = YEAR, y = PREVALENCE, group = as.factor(CATEGORY))) +
            geom_point(aes(colour = CATEGORY), size = 3) +
            geom_line(aes(colour = CATEGORY), size = 1) +
            scale_x_continuous(limits = c(1988, current_year),
                               breaks = seq(1990, current_year, 5)) +
            scale_color_manual(labels = c("Immature prevalence", "Mature prevalence", 
                                          "Population prevalence"), 
                               values = c("#E69F00", "#56B4E9", "#009E73")) +
            labs(y = "Disease Prevalence (%)", x = "") +
            theme_bw() +
            theme(legend.text = element_text(size = 11),
                  axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 14), 
                  axis.text.y = element_text(size = 12),
                  legend.title = element_blank())


# Faceted plot 
labels <- function(variable,value){
  return(names[value])
}

names <- list("IMM_PREVALENCE" = "Immature prevalence", 
              "MAT_PREVALENCE" = "Mature prevalence", 
              "POP_PREVALENCE" = "Population prevalence")
facet_plot <- ggplot(plot_dat, aes(x = YEAR, y = PREVALENCE, group = as.factor(CATEGORY))) +
              geom_point(aes(colour = CATEGORY), size = 3) +
              geom_line(aes(colour = CATEGORY), size = 1) +
              scale_x_continuous(limits = c(1988, current_year),
                                 breaks = seq(1990, current_year, 5)) +
              scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
              labs(y = "Disease Prevalence (%)", x = "") +
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
imm_plot <- ggplot(plot_dat %>% filter(CATEGORY == "IMM_PREVALENCE"), 
                   aes(x = YEAR, y = PREVALENCE)) +
            geom_point(colour = "#009E73", size = 3) +
            geom_line(colour = "#009E73", size = 1) +
            scale_x_continuous(limits = c(1988, current_year),
                               breaks = seq(1990, current_year, 5)) +
            geom_hline(aes(yintercept = mean(PREVALENCE, na.rm = TRUE)), linetype = 5) +
            labs(y = "Immature Disease Prevalence (%)", x = "") +
            theme_bw() +
            theme(axis.title.y = element_text(size = 14),
                  axis.text.x = element_text(size = 12), 
                  axis.text.y = element_text(size = 12),
                  legend.title = element_blank())
  
ggsave("./figures/bcd_imm_prev.png", imm_plot,
       height = 6, width = 10)

