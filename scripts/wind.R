## WIND
## 60° component of the average May - June wind vector (need in m/s?)
## - back-calculate u and v components of hourly vectors
## - calculate 60 deg vector component by rotating angle??




## NCEP - Tyler Hennon
wind <- read.csv(paste0(data_dir, "wind_NCEP.csv")) %>% 
        rename(ALONG_SHELF = ALONG.SHELF.WINDS..M.S.,
               CROSS_SHELF = CROSS.SHELF.WINDS..M.S.,
               year = YEAR) %>%
        filter(MONTH %in% c(5:6)) %>%
        group_by(year) %>%
        summarise(wind_along_shelf = mean(ALONG_SHELF, na.rm = TRUE),
                  wind_cross_shelf = mean(CROSS_SHELF, na.rm = TRUE)) %>%
        write_csv("./outputs/wind_NCEP.csv")


# Plot
ggplot(wind, aes(x = year, y = wind_along_shelf)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Mean May-June along-shelf wind (m/s)") +
  geom_hline(aes(yintercept = mean(wind_along_shelf, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(wind_along_shelf, na.rm = TRUE) - sd(wind_along_shelf, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(wind_along_shelf, na.rm = TRUE) + sd(wind_along_shelf, na.rm = TRUE)), color = "green4") +
  # xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "wind_along_shelf.png"), height =  4, width = 6)

ggplot(wind, aes(x = year, y = wind_cross_shelf)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Mean May-June cross-shelf wind (m/s)") +
  geom_hline(aes(yintercept = mean(wind_cross_shelf, na.rm = TRUE)), linetype = 5) +
  geom_hline(aes(yintercept = mean(wind_cross_shelf, na.rm = TRUE) - sd(wind_cross_shelf, na.rm = TRUE)), color = "green4") +
  geom_hline(aes(yintercept = mean(wind_cross_shelf, na.rm = TRUE) + sd(wind_cross_shelf, na.rm = TRUE)), color = "green4") +
  # xlim(min(years), max(years)) +
  theme_bw()
ggsave(paste0(fig_dir, "wind_cross_shelf.png"), height =  4, width = 6)



# ## ASOS ------------
# wind <- read.csv(paste0(data_dir, "wind_ASOS.csv")) %>% 
#         mutate(date = as.Date(valid),
#                adj_dir = drct - 60, 
#                u = sknt*cos(adj_dir)/0.51444) %>%
#         separate(date, into = c("year", "month", "day")) %>%
#         filter(month %in% c("05", "06")) %>%
#         mutate(year = as.numeric(year)) %>%
#         group_by(day, month, year) %>%
#         summarise(mean_60deg = mean(u, na.rm = TRUE)) %>%
#         # ungroup() %>%
#         group_by(month, year) %>%
#         summarise(mean_60deg = mean(mean_60deg, na.rm = TRUE)) %>%
#         group_by(year) %>%
#         summarise(mean_60deg = mean(mean_60deg, na.rm = TRUE))
#   
# 
# 
# # Plot
# ggplot(wind, aes(x = year, y = mean_60deg)) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Year", y = "Mean summer 60 degree wind component (m/s)") +
#   geom_hline(aes(yintercept = mean(mean_60deg, na.rm = TRUE)), linetype = 5) +
#   # xlim(min(years), max(years)) +
#   theme_bw()
# ggsave(paste0(fig_dir, "wind.png"), height =  4, width = 6)
# 
# 
# # # Save output
# # wind %>%
# #   write_csv("./outputs/wind.csv")