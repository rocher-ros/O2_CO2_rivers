# Information ----
# Script to assess patterns of CO2-O2 in river versus lakes
# Author: Gerard Rocher-Ros
# Last edit: 2024-10-18

# Install and Load libraries ----

# List of all packages needed
package_list <- c('tidyverse', 'patchwork', 'readxl')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if (length(packages_missing) >= 1) {
  install.packages(packages_missing)
}

# Now load all the packages
lapply(package_list, require, character.only = TRUE)


# Functions ----
# Function to estimate ellipse properties, from Vachon et al . 2020
elipse <- function(x, y) {
  if (length(x) > 50) {
    corMat = cor(cbind(x, y))
    covMat <- var(cbind(x, y))
    evals <- eigen(covMat)$values
    ell.len <- 2 * sqrt(5.991 * evals)

    out <- tibble(width = ell.len[2], length = ell.len[1])
  } else {
    out <- tibble(width = NA, length = NA)
  }
  out
}

# Read files ----
## read river data  ----

#site data from streampulse,
sp_site_labels <- read_csv(
  "prepared data/river data/streampulse/streampulse_site_data.csv"
) %>%
  filter(str_detect(variableList, "CO2_ppm"))

sp_site_data <- read_csv(
  "prepared data/river data/streampulse/site_info_streams.csv"
) %>%
  select(-geometry)

#metab data from streampulse
sp_metab_data <- read_csv(
  "prepared data/river data/streampulse/all_daily_model_results.csv"
) %>%
  filter(GPP >= 0, ER <= 0, valid_day == 1) %>%
  mutate(day_month = as.Date(date))

sp_metab_avg <- sp_metab_data %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .by = "site")

#discharge data from USGS of all sites required, join with site info
discharge_all <- read_csv(
  "prepared data/river data/streampulse/discharge_usgs.csv"
) %>%
  mutate(date = as.Date(date_time)) %>%
  summarise(discharge = mean(discharge), .by = c("site", "date")) %>%
  left_join(
    sp_site_labels %>% select(siteID, USGSgageID),
    by = c("site" = "USGSgageID")
  ) %>%
  mutate(discharge_m3s = discharge * 0.028316832)

#read co2_o2 data, preprocessed from before
sp_data <- read_csv(
  "prepared data/river data/streampulse/stream_dataset.csv"
) %>%
  mutate(date = as.Date(dateTimeUTC)) %>%
  left_join(discharge_all, by = c("siteID", "date")) %>%
  mutate(
    Discharge_joined = ifelse(
      is.na(Discharge_m3s) == TRUE,
      discharge_m3s,
      Discharge_m3s
    ),
    siteID = str_remove(siteID, "/")
  ) %>%
  left_join(sp_metab_data, by = c("siteID" = "site", "date" = "day_month"))


river_coords <- read_excel("prepared data/river data/table_s1.xlsx") %>%
  mutate(site = `Site ID`, type = "Rivers") %>%
  select(type, site, Latitude, Longitude)


## read lake data ----
# read the data deposited in Vachon et al. 2020

lakes_coords <- read_csv(
  "prepared data/lake data/vachon2020/table1_sensor_info.csv"
) %>%
  mutate(type = "Lakes") %>%
  select(type, site = Lake, Latitude, Longitude)


files_lakes <- list.files(
  "prepared data/lake data/vachon2020",
  pattern = "*.csv",
  full.names = T
)

lakes <- files_lakes %>%
  set_names() %>%
  map_df(
    ~ read_delim(.x, delim = ";", escape_double = FALSE, trim_ws = TRUE) %>%
      mutate(file_name = basename(.x))
  ) %>%
  mutate(site = file_name %>% str_remove(".csv")) %>%
  select(-file_name)


unique(lakes$site)


sites_coords <- bind_rows(lakes_coords, river_coords)

# Data exploration and pre-processing ----

#lakes
lakes %>%
  ggplot(aes(CO2dep, O2dep)) +
  geom_density_2d() +
  geom_abline(slope = -1, intercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)

#rivers
sp_data %>%
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3)) +
  stat_ellipse(aes(color = siteID)) +
  geom_abline(slope = -1, intercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)

#Join both datasets
rivers_lakes <- lakes %>%
  mutate(type = "Lakes") %>%
  select(type, CO2dep, O2dep, site) %>%
  bind_rows(
    sp_data %>%
      mutate(type = "Rivers", CO2dep = CO2dep_mmolm3, O2dep = O2dep_mmolm3) %>%
      select(type, CO2dep, O2dep, site = siteID)
  ) %>%
  drop_na(CO2dep, O2dep) %>%
  left_join(sites_coords)

#summary stats
rivers_lakes %>% group_by(type) %>% count(site) %>% print(n = 50)
#12 lakes and 30 streams

# Main plots ----
densities_2dplot <-
  rivers_lakes %>%
  ggplot(aes(CO2dep, O2dep)) +
  geom_abline(slope = -1, intercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_hline(yintercept = 0, color = "gray40") +
  stat_ellipse(
    aes(color = type),
    geom = "polygon",
    alpha = 0,
    level = 0.9,
    linewidth = .7
  ) +
  stat_density_2d_filled(
    aes(fill = type, color = type, group = type, alpha = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    contour_var = "ndensity",
    breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1),
    linewidth = .2
  ) +
  #geom_point(data=. %>% summarise(CO2dep = median(CO2dep), #in case you want to see all the points, takes some time
  #                                O2dep = median(O2dep), .by= type),
  #           aes(fill=type), shape= 23, size= 3, show.legend = FALSE)+
  scale_y_continuous(breaks = c(-150, -100, -50, 0, 50), limits = c(-180, 80)) +
  scale_x_continuous(
    breaks = c(-50, 0, 50, 100, 150, 200, 250),
    limits = c(-60, 250)
  ) +
  scale_fill_manual(values = c("dodgerblue3", "darkorange")) +
  scale_color_manual(values = c("dodgerblue3", "darkorange")) +
  labs(
    x = expression(CO[2] ~ departure ~ (mu * mol ~ L^-1)),
    y = expression(O[2] ~ departure ~ (mu * mol ~ L^-1)),
    color = "",
    fill = ""
  ) +
  guides(alpha = "none") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    panel.border = element_rect(linewidth = 1, fill = NA),
    plot.margin = margin(-20, -20, 0, 0)
  )


co2_density_plot <- rivers_lakes %>%
  ggplot() +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_density(
    aes(CO2dep, y = after_stat(scaled), fill = type, color = type),
    alpha = .3
  ) +
  geom_point(
    data = . %>% summarise(CO2dep = median(CO2dep), .by = type),
    aes(x = CO2dep, y = 0.2, color = type),
    shape = 19,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("dodgerblue3", "darkorange")) +
  scale_color_manual(values = c("dodgerblue3", "darkorange")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(
    breaks = c(-50, 0, 50, 100, 150, 200, 250),
    limits = c(-60, 250)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(5, 5, -20, 5)
  ) +
  labs(x = "", y = "")

o2_density_plot <- rivers_lakes %>%
  ggplot() +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_density(
    aes(O2dep, y = after_stat(scaled), fill = type, color = type),
    alpha = .3
  ) +
  geom_point(
    data = . %>% summarise(O2dep = median(O2dep), .by = type),
    aes(x = O2dep, y = 0.2, color = type),
    shape = 19,
    size = 3,
    show.legend = FALSE
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("dodgerblue3", "darkorange")) +
  scale_color_manual(values = c("dodgerblue3", "darkorange")) +
  scale_x_continuous(breaks = c(-150, -100, -50, 0, 50), limits = c(-180, 80)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(5, 5, -20, 5)
  ) +
  coord_flip() +
  labs(x = "", y = "")


plot_lakes_rivers <- co2_density_plot +
  guide_area() +
  densities_2dplot +
  o2_density_plot +
  plot_layout(
    ncol = 2,
    nrow = 2,
    guides = 'collect',
    widths = c(5, 1),
    heights = c(1, 4)
  )

plot_lakes_rivers

ggsave(
  filename = "plots/main/fig1_lakes_rivers.png",
  plot_lakes_rivers,
  scale = 1,
  height = 6,
  width = 6.5
)


# Calculate ellipse metrics ----
#this is used for the overview in the introduction
metrics_site <- rivers_lakes %>%
  drop_na(CO2dep, O2dep) %>%
  group_by(type) %>%
  summarise(
    meanCO2 = mean(CO2dep),
    meanO2 = mean(O2dep),
    n = n(),
    offset = meanCO2 + meanO2,
    elipse = elipse(CO2dep, O2dep)
  ) %>%
  unnest(elipse)

metrics_site
