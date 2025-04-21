## WIND
## 60° component of the average May - June wind vector (need in m/s?)
## - back-calculate u and v components of hourly vectors
## - calculate 60 deg vector component by rotating angle??

wind <- read.csv(paste0(data_dir, "wind_ASOS.csv")) %>% 
        mutate(date = as.Date(valid),
               adj_dir = drct - 60, 
               u = sknt*cos(adj_dir)/0.51444) %>%
        separate(date, into = c("year", "month", "day")) %>%
        filter(month %in% c("05", "06")) %>%
        mutate(year = as.numeric(year)) %>%
        group_by(day, month, year) %>%
        summarise(mean_60deg = mean(u, na.rm = TRUE)) %>%
        # ungroup() %>%
        group_by(month, year) %>%
        summarise(mean_60deg = mean(mean_60deg, na.rm = TRUE)) %>%
        group_by(year) %>%
        summarise(mean_60deg = mean(mean_60deg, na.rm = TRUE))
  


# Plot
ggplot(wind, aes(x = year, y = mean_60deg)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Mean summer 60 degree wind component (m/s)") +
  geom_hline(aes(yintercept = mean(mean_60deg, na.rm = TRUE)), linetype = 5) +
  # xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "wind.png"), height =  4, width = 6)


# Save output
mean_size %>%
  rename(year = YEAR,
         female_sam = MEAN_SIZE) %>%
  write_csv("./outputs/female_SAM.csv")