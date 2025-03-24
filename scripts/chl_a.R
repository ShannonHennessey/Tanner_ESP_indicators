## Purpose: To calculate April - July mean OCCCI chlorophyll-a concentration in
##          BISERP regions 3-6,8
##
## NOTES:
## - think about which time period, which regions
##    - larvae in water column through July, but maybe earlier is more significant?
##    - mat female distribution (3-5) vs. potential larval advection (6,8 too)?
## - eventually add in peak bloom timing -- some sort of DFA/combined indicator?


## Read in setup
source("./scripts/setup.R")


## Get extent of BSIERP regions
bsierp <- akgfmaps::get_bsierp_regions(set.crs = 4326)

ext <- bsierp %>%
       filter(BSIERP_ID %in% c(3,4,5,6,8)) %>%
       st_bbox()
 


# ## Download data
# for(i in years){
#   file_name <- paste0("./data/nc/occ8_", i ,".nc")
#   download.file(url = paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/pmlEsaCCI60OceanColor8Day.nc?chlor_a%5B(",
#                              i,"-04-01T00:00:00Z):(", i,"-07-31T00:00:00Z)%5D%5B(",ext[2],"):(",ext[4],")%5D%5B(",ext[1],
#                              "):(",ext[3],")%5D&.draw=surface&.vars=longitude%7Clatitude%7Cchlor_a&.colorBar=%7C%7C%7C%7C%7C&.bgColor=0xffccccff"),
#                 method = "libcurl", mode = "wb", destfile = file_name)
# }


## Load/compile data
tidy_chl <- function(file) {
  tidync(file) %>% 
    hyper_tibble() %>% 
    mutate(date = as_datetime(time),
           latc = as.numeric(substr(latitude, 1, 7)),
           longc = as.numeric(substr(longitude, 1, 9)),
           chlorophyll = round(chlor_a, 3),
           year = year(date),
           month = month(date),
           day = day(date))
}

# Append all files
datalist = list()

for(i in years){
  dat <- tidy_chl(paste0("./data/nc/occ8_", i ,".nc"))
  datalist[[i]] <- dat 
}

# Convert list to data frame
occ <- dplyr::bind_rows(datalist)

# from https://github.com/MattCallahan-NOAA/chla-indicator-comparison/blob/main/chla-indicator-comparison-data.RMD
occ_grid <- readRDS("./data/occ_chl_spatial_lookup.RDS")
occ_esp <- occ %>%
           inner_join(occ_grid, by = c("longc" = "longitude", "latc" = "latitude"))


# Plot BSIERP regions
occ_esp %>%
  filter(BSIERP_ID %in% c(3, 4, 5, 6, 8),
         month %in% c(4:7)) %>%
  group_by(year) %>%
  summarize(mean_chla = mean(chlorophyll, na.rm = T)) %>%
  ggplot() +
  geom_line(aes(x = year, y = mean_chla)) +
  geom_point(aes(x = year, y = mean_chla)) +
  labs(x = "Year", y = "Mean Chlorophyll-a Concentration (µg/L)") +
  xlim(1988, 2024) +
  theme_bw() 
ggsave("./figures/BSIERP_OCCCI_AMJJ.png")






