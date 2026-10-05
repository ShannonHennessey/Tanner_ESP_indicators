## Purpose: To visualize Tanner densities/distributions and cold pool extent
##
## Authors: S. Hennessey, E. Fedewa, E. Ryznar
##
## NOTES: 
##


## Load packages
library(tidyverse)
# library(tidync)
# library(sf)
# library(terra)
library(akgfmaps)
# library(ggridges)
# library(patchwork)
# library(hrbrthemes)
# library(ggtext)
# library(ggpubr)
# library(ggh4x)
# library(gridExtra)
# library(lemon)

#install.packages("remotes")
#remotes::install_github("afsc-gap-products/akgfmaps")

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

## Extract haul/temperature data
haul <- tanner$haul


## Load spatial data
in.crs <- "+proj=longlat +datum=NAD83" # CRS is in lat/lon
map.crs <- "EPSG:3338" # final crs for mapping/plotting: Alaska Albers
my_colors <- c("#D55E00","#a6bddb", "#74a9cf", "#0570b0", "#034e7b") # plot colors

survey_strata <- terra::vect(paste0(data_dir, "SAP_layers.gdb"), layer = "EBS_grid")
# #EBS/NBS Boundary line
# boundary <- st_read(layer = "EBS_NBS_divide", paste0(data_dir, "SAP_layers.gdb")) 

ebs_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")
ebs_survey_areas <- ebs_layers$survey.area
ebs_survey_areas$survey_name <- c("Eastern Bering Sea", "Northern Bering Sea")

# #Survey areas plot
# ggplot() +
#   geom_sf(data = ebs_layers$akland) +
#   geom_sf(data = ebs_survey_areas, mapping = aes(fill = survey_name)) +
#   scale_x_continuous(limits = ebs_layers$plot.boundary$x,
#                      breaks = ebs_layers$lon.breaks) +
#   scale_y_continuous(limits = ebs_layers$plot.boundary$y,
#                      breaks = ebs_layers$lat.breaks) +
#   scale_fill_viridis_d(name = "Survey") +
#   theme_bw()
# 
# #EBS/NBS shelf Survey grid plot
# ggplot() +
#   geom_sf(data = ebs_layers$akland) +
#   geom_sf(data = ebs_layers$survey.grid, fill = NA) +
#   # geom_sf(data= boundary, linewidth = 2) +
#   scale_x_continuous(limits = ebs_layers$plot.boundary$x,
#                      breaks = ebs_layers$lon.breaks) +
#   scale_y_continuous(limits = ebs_layers$plot.boundary$y,
#                      breaks = ebs_layers$lat.breaks) +
#   theme_bw()




## Compute mean summer bottom temperature
bottom_temp <- haul %>%
               filter(YEAR %in% years,
                      !HAUL_TYPE == 17) %>%
               distinct(YEAR, STATION_ID, GEAR_TEMPERATURE, MID_LATITUDE, MID_LONGITUDE) %>%
               # group_by(YEAR) %>%
               # summarise(summer_bt = mean(GEAR_TEMPERATURE, na.rm = T)) %>%
               # right_join(expand.grid(YEAR = years)) %>%
               mutate(COLDPOOL = ifelse(GEAR_TEMPERATURE < 2, 1, 0)) %>%
               # transform into spatial data frame
               sf::st_as_sf(coords = c("MID_LONGITUDE", "MID_LATITUDE"), crs = in.crs) %>%
               sf::st_transform(sf::st_crs(map.crs)) %>%
               vect(.) %>%
               # mask(., survey_strata) %>% # limit to survey region
               sf::st_as_sf() 



# Transform crab data into spatial data frame
crab_dat <- cpue %>% 
            # Convert lat/long to an sf object
            st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(4326)) %>%
            #st_as_sf needs crs of the original coordinates- need to transform to Alaska Albers
            st_transform(crs = map.crs)

# Use the spatial data frame to generate a convex hull around the data extent
coldpool_hull <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(bottom_temp %>% filter(COLDPOOL == 1)))), 
                                  dist = 15000), dTolerance = 5000)


#Add hull, ice extent layer and crab data to map
base_map <- ggplot() +
            geom_sf(data = ebs_layers$survey.grid, fill = NA, color = alpha("grey80")) +
            geom_sf(data = ebs_survey_areas, fill = NA) +
            geom_sf(data = ebs_layers$akland, fill = "grey80", color = "black") +
            #add hull for sea ice extent- though we won't use this approach here....
            #geom_sf(data = ice_hull,
            #fill = NA,
            #color = alpha("red", 0.85),
            #linewidth = 1) +
            # add bottom temperature
            geom_sf(data = bottom_temp, mapping = aes(color = GEAR_TEMPERATURE), alpha = 0.25) +
            # add crab
            # geom_sf(crab_dat, mapping = aes(size = CPUE), alpha = 0.6) +
            scale_x_continuous(limits = ebs_layers$plot.boundary$x,
                               breaks = ebs_layers$lon.breaks) +
            scale_y_continuous(limits = ebs_layers$plot.boundary$y,
                               breaks = ebs_layers$lat.breaks) +
            scale_size_continuous(range = c(1, 4)) +
            theme_bw() +
            facet_wrap(~YEAR) +
            scale_color_manual(values = c("#034e7b", "#238b45")) +
            theme(legend.position = "bottom",
                  legend.margin = margin(-5, 0, -1, 0), # reducing white space b/w plot and legend
                  legend.spacing.x = unit(-2, "mm"),
                  legend.spacing.y = unit(-2, "mm")) +
            theme(plot.margin = margin(0, -5, 0, -5)) +
            theme(axis.text = element_text(size = 8)) +
            theme(axis.text.x = element_blank())

#Goofy workaround for duplicating legends by lme:
#Add EBS-specific legend and extract 
ebs_map <- base_map +
           labs(x = "", y = "", size = "CPUE") +
           guides(color = "none") +
           guides(size = guide_legend(override.aes = list(color="#034e7b"))) +
           theme(legend.title = element_text(size = 9))
#extract EBS legend
ebs_legend <- g_legend(ebs_map)

#Add NBS-specific legend and extract 
base_map +
  labs(x="", y="", size = expression(paste("Northern Bering Sea samples"))) +
  guides(color = "none") +
  guides(size = guide_legend(override.aes = list(color="#238b45"))) +
  theme(legend.title=element_text(size=9)) -> nbs_map
#extract NBS legend
nbs_legend <- g_legend(nbs_map)

#Base map with workaround to center bottom panel 
base_map +
  guides(color = "none", size = "none") -> base
#workaround to center the bottom panel
design <- c(
  "
AABBCC
#DDEE#
"
)
base + ggh4x::facet_manual(~year, design=design) -> final_base

#Now combine map and two legends for final map
final_base / nbs_legend / ebs_legend +
  plot_layout(heights= c(5,.1,.1)) -> final_map

### PANEL A ------------------------------------------------------------------
#EBS and NBS abundance timeseries 

#calculate EBS abundance timeseries 
ebs_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(HAUL_TYPE == 3,
         YEAR >= 1988) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue_cnt = ncrab / AREA_SWEPT) %>%
  # join to hauls that didn't catch crab 
  right_join(ebs_haul %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(HAUL_TYPE ==3,
                      YEAR >= 1988) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue_cnt = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  
  #join to stratum
  left_join(ebs_strata %>%
              select(STATION_ID, SURVEY_YEAR, STRATUM, TOTAL_AREA) %>%
              filter(SURVEY_YEAR >= 1988) %>%
              rename_all(~c("GIS_STATION", "YEAR",
                            "STRATUM", "TOTAL_AREA"))) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA) %>%
  summarise(MEAN_CPUE = mean(cpue_cnt , na.rm = T),
            N_CPUE = n(),
            VAR_CPUE = (var(cpue_cnt)*(TOTAL_AREA^2))/N_CPUE,
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
  distinct() %>%
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6,
            SD_CPUE = sqrt(sum(VAR_CPUE)),
            ABUNDANCE_CI = (1.96*(SD_CPUE))/1e6) %>%
  bind_rows(missing <- data.frame(YEAR = 2020)) %>%
  mutate(lme = rep("Collapsing Eastern Bering Sea")) %>%
  arrange(YEAR) -> ebs_abundance

#----------------------
#calculate NBS abundance timeseries 
nbs_haul %>%
  mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
  filter(YEAR >= 1988) %>%
  group_by(YEAR, GIS_STATION, AREA_SWEPT) %>%
  summarise(ncrab = sum(SAMPLING_FACTOR, na.rm = T)) %>%
  ungroup %>%
  # compute cpue per nmi2
  mutate(cpue_cnt = ncrab / AREA_SWEPT) %>%
  # join to hauls that didn't catch crab 
  right_join(nbs_haul %>% 
               mutate(YEAR = as.numeric(str_extract(CRUISE, "\\d{4}"))) %>%
               filter(YEAR >= 1988) %>%
               distinct(YEAR, GIS_STATION, AREA_SWEPT)) %>%
  replace_na(list(cpue_cnt = 0)) %>%
  replace_na(list(ncrab = 0)) %>%
  
  #join to stratum
  left_join(nbs_strata %>%
              select(GIS_STATION, SURVEY_YEAR, STRATUM, TOTAL_AREA) %>%
              filter(SURVEY_YEAR >= 1988) %>%
              rename_all(~c("GIS_STATION", "YEAR",
                            "STRATUM", "TOTAL_AREA"))) %>%
  #Scale to abundance by strata
  group_by(YEAR, STRATUM, TOTAL_AREA) %>%
  summarise(MEAN_CPUE = mean(cpue_cnt , na.rm = T),
            N_CPUE = n(),
            VAR_CPUE = (var(cpue_cnt)*(TOTAL_AREA^2))/N_CPUE,
            ABUNDANCE = (MEAN_CPUE * mean(TOTAL_AREA))) %>%
  distinct() %>%
  filter(STRATUM != "NA") %>% #NAs are issues documented in 11/17 SAP email! 
  group_by(YEAR) %>%
  #Sum across strata
  summarise(ABUNDANCE_MIL = sum(ABUNDANCE)/1e6,
            SD_CPUE = sqrt(sum(VAR_CPUE)),
            ABUNDANCE_CI = (1.96*(SD_CPUE))/1e6) %>%
  mutate(lme = rep("Non-collapsing Northern Bering Sea")) -> nbs_abundance

#Combine EBS and NBS datasets and plot ------------------------------------------
ebs_abundance %>%
  bind_rows(nbs_abundance) -> plot_abun

#Plot with lines
plot_abun %>%
  filter (YEAR >= 1995) %>%
  ggplot() +
  geom_point(data = subset(plot_abun, lme == "Collapsing Eastern Bering Sea"), 
             aes(y=ABUNDANCE_MIL, x=YEAR, color = lme), size = 2) +
  geom_line(data = subset(plot_abun, lme == "Collapsing Eastern Bering Sea"), 
            aes(y=ABUNDANCE_MIL, x=YEAR, color = lme), size = 1, alpha = 7) +
  geom_ribbon(data = subset(plot_abun, lme == "Collapsing Eastern Bering Sea"),
              aes(ymin = ABUNDANCE_MIL - ABUNDANCE_CI, 
                  ymax = ABUNDANCE_MIL + ABUNDANCE_CI, x=YEAR),
              fill = "#9ECAE1", alpha = 0.2) +
  geom_point(data = subset(plot_abun, lme == "Non-collapsing Northern Bering Sea"), 
             aes(y=ABUNDANCE_MIL, x=YEAR, color = lme), size = 2) +
  geom_errorbar(data = subset(plot_abun, lme == "Non-collapsing Northern Bering Sea"),
                aes(ymin = ABUNDANCE_MIL - ABUNDANCE_CI, 
                    ymax = ABUNDANCE_MIL + ABUNDANCE_CI, x=YEAR),
                color = "#74C476", alpha = .9) +
  scale_color_manual(values = c("#4292C6", "#41AB5D")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(y="Snow Crab \nAbundance (millions)", x="") +
  theme(legend.position=c(.2,.13)) +
  theme(legend.background = element_rect(color="transparent", fill="transparent")) +
  theme(axis.title.y = element_text(size=10)) +
  theme(legend.title=element_blank()) +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(limits=c(350,NA), breaks=c(500,1500, 5000,10000,20000)) +
  theme(legend.text=element_text(size=9.5)) +
  guides(color = guide_legend(override.aes = list(size=2.5, shape=19))) +
  theme(plot.margin = margin(0,0.5,0,0)) -> abun_plot

#Stacked bar plot
plot_abun %>%
  filter (YEAR >= 1992) %>%
  ggplot(aes(fill=lme, y=ABUNDANCE_MIL, x=YEAR)) + 
  geom_col(position = position_dodge(preserve = 'single')) +
  theme_bw() +
  scale_fill_manual(values = c("#4292C6", "#084594")) +
  labs(y="Snow Crab Abundance (millions)", x="") +
  theme(legend.title=element_blank())

### COMBINE PANELS AND SAVE FIGURE --------------------------------------------

#Figure 1 for ms: combined abun, and map plot
abun_plot + plot_annotation(tag_levels = 'a') 
ggsave(paste0(fig_dir, "Fig1a.png"), height=4 , width=7.5, units="in")

final_map + plot_annotation(tag_levels = list('b')) &
  theme(plot.tag.position  = c(.02, 1))
ggsave(paste0(fig_dir, "Fig1b.png"), height=6 , width=7, units="in")
#These were manually combined in pwpt as patchwork was distorting map size! 

### FIG 2 -----------------------------------------------------------
#Barplots of annual mean cpue and temp for EBS and NBS

#cpue data wrangling
condition_master %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  mutate(lme = recode(lme, EBS = "Collapsing Eastern Bering Sea", NBS = "Non-collapsing Northern Bering Sea")) %>%
  distinct(year, lme, gis_station, cpue) %>%
  group_by(year, lme) %>%
  summarize(mean_cpue = mean(cpue, na.rm = T)/1000, #converting to thous crab/nmi2
            sd_cpue = sd(cpue, na.rm = T)/1000,
            n_cpue = n()) %>%
  mutate(se_cpue = sd_cpue / sqrt(n_cpue),
         lower.ci = mean_cpue - qt(1 - (0.05 / 2), n_cpue - 1) * se_cpue,
         upper.ci = mean_cpue + qt(1 - (0.05 / 2), n_cpue - 1) * se_cpue,
         year = as.factor(year)) -> cpue.dat
#and plot
ggplot(cpue.dat, aes(year, mean_cpue)) +
  geom_col(aes(fill = ordered(year)), size=3) +
  geom_errorbar(aes(year, ymin=mean_cpue - se_cpue, ymax=mean_cpue + se_cpue), 
                width=0.3, size=0.5, color = "grey40") +
  labs(y = "Mean Snow Crab Density<br>(thousand crab/nmi<sup> 2</sup>)\n", x= "") +
  scale_fill_manual(values=my_colors) +
  facet_wrap(~lme, labeller = label_wrap_gen(multi_line = TRUE)) +
  geom_vline(data = subset(cpue.dat, lme == "Collapsing Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(cpue.dat, lme == "Collapsing Eastern Bering Sea"), aes(x = 1, y=700, label = "Mid-collapse"),
            size = 2.2, color = "#D55E00") +
  geom_text(data = subset(cpue.dat, lme == "Collapsing Eastern Bering Sea"), aes(x = 3, y=700, label = "Post-collapse"),
            size = 2.2, color = "#0072B2") +
  theme_ipsum(axis_title_size = 10.5, axis_text_size =10) +
  theme(legend.position="none") +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(limits=c(0,NA), breaks=c(5, 25,75,200,600)) +
  theme(axis.title.y=element_text(colour="grey30", size = 10.5)) +
  theme(strip.text = element_text(hjust = .5, size=11)) +
  theme(axis.title.y = element_textbox_simple(orientation = "left-rotated", halign = 0.5)) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())  -> mean_cpue_plot

#temperature data wrangling
condition_master %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  mutate(lme = recode(lme, EBS = "Collapsing Eastern Bering Sea", NBS = "Non-collapsing Northern Bering Sea")) %>%
  distinct(year, lme, gis_station, gear_temperature) %>%
  group_by(year, lme) %>%
  summarize(mean_temp = mean(gear_temperature, na.rm = T),
            sd_temp = sd(gear_temperature, na.rm = T),
            n_temp = n()) %>%
  mutate(se_temp = sd_temp / sqrt(n_temp),
         lower.ci = mean_temp - qt(1 - (0.05 / 2), n_temp - 1) * se_temp,
         upper.ci = mean_temp + qt(1 - (0.05 / 2), n_temp - 1) * se_temp,
         year = as.factor(year)) -> temp.dat
#plot  
ggplot(temp.dat, aes(year, mean_temp)) +
  geom_col(aes(fill = ordered(year)), size=3) +
  geom_errorbar(aes(year, ymin = ifelse(mean_temp - se_temp < 0, 0, mean_temp - se_temp), 
                    ymax=mean_temp + se_temp),
                width=0.3, size=0.5, color = "grey40") +
  labs(y = expression("Mean Bottom \nTemperature" ( degree~C)), x = "") +
  scale_fill_manual(values=my_colors) +
  facet_wrap(~lme) +
  geom_vline(data = subset(cpue.dat, lme == "Collapsing Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(cpue.dat, lme == "Collapsing Eastern Bering Sea"), aes(x = 1, y=3.3, label = "Mid-collapse"),
            size = 2.2, color = "#D55E00") +
  geom_text(data = subset(cpue.dat, lme == "Collapsing Eastern Bering Sea"), aes(x = 3, y=3.3, label = "Post-collapse"),
            size = 2.2, color = "#0072B2") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 10.5, axis_text_size =10) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(colour="grey30", hjust = 0.5, size = 10.5)) +
  theme(strip.text = element_blank()) +
  theme(panel.grid.major.x = element_blank())  -> mean_temp_plot

### COMBINE PANELS AND SAVE FIGURE --------------------------------------------

#Figure 2 for ms: combined density and temperature plot
mean_cpue_plot / plot_spacer() / mean_temp_plot  + plot_layout(heights = c(6, -3 , 6)) +
  plot_annotation(tag_levels = list(c('a', 'b'))) 
ggsave(paste0(fig_dir, "Fig2.png"), height=7 , width=7.5, units="in")


-----------------------------------------------------------------------------
  #Bonus Figs: And now pdfs by year
  
  #temperature
  condition_master  %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  group_by(year, lme, gis_station) %>%
  summarise(temperature = mean(gear_temperature)) %>%
  ggplot(aes(temperature,factor(year))) +
  geom_density_ridges(aes(fill=factor(year)), scale=2,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x),
                      rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2,
                      alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "Bottom Temperature (C)", y = "Count") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) 

#cpue 
condition_master  %>%
  filter(lme %in% c("EBS", "NBS")) %>%
  group_by(year, lme, gis_station) %>%
  summarise(cpue = mean(cpue)) %>%
  ggplot(aes(cpue,factor(year))) +
  geom_density_ridges(aes(fill=factor(year)), scale=2,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x),
                      rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2,
                      alpha = 0.7) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~lme) +
  theme_bw() +
  labs(x= "Snow Crab Density", y = "Count") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) 


condition_master %>%
  group_by(year, lme) %>%
  summarise(mean_cpue = mean(cpue), 
            mean_temp = mean(gear_temperature)) %>%
  ggplot( ) +
  geom_point(aes(mean_cpue, mean_temp, color = lme)) +
  geom_line(aes(mean_cpue, mean_temp, color = lme)) +
  geom_text(label = year)

condition_master %>%
  mutate(fourth.root.cpue = as.numeric(cpue^0.25)) %>%
  ggplot( ) +
  geom_point(aes(fourth.root.cpue, gear_temperature, color = lme)) +
  geom_line(aes(fourth.root.cpue, gear_temperature, color = lme)) +
  geom_text(label = year)