# Load  and install packages ----
# List of all packages needed
package_list <- c('tidyverse', 'ellipse', 'lmodel2')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

# Now load all the packages
lapply(package_list, require, character.only = TRUE)


theme_set(theme_bw())

#functions to calculate ellipse and regression metrics ----
# Modified from Vachon et al., 2020 GCB

elipse <- function(x, y){
  
  if(length(x) > 50){
    corMat = cor(cbind(x,y))
    covMat <- var(cbind(x,y))
    evals <- eigen(covMat)$values
    ell.len <- 2*sqrt(5.991*evals)
    
    out <- tibble(width = ell.len[2],
                  length = ell.len[1]) 
  }else{
    out <- tibble(width = NA,
                  length = NA) }
  out
}

#Type II regression for slope calculation 
regression_metrics <- function(x, y){
  
  if(length(x) > 50){  
    reg <- lmodel2(y ~ x, nperm = 99)
    slope <- reg$regression.results[2,3]
    corr <- reg$r
    
    out <- tibble(EQ= 1/abs(slope), corr= corr) 
  }else{
    out <- tibble(EQ = NA,
                  corr = NA) }
  out
  
}

# Read files ----
sp_site_labels <- read_csv("data/stream_datasets/_raw/streampulse/streampulse_site_data.csv") %>% 
  filter(str_detect(variableList,"CO2_ppm" ))

sp_site_data <- read_csv("data/stream_datasets/site_info_streams.csv") %>% 
  select(-geometry)


sp_metab_data <- read_csv("data/stream_datasets/all_daily_model_results.csv") %>% 
  filter(GPP>= 0, ER <= 0, valid_day == 1) %>% 
  mutate(day_month= as.Date(date))

sp_metab_avg <- sp_metab_data %>% 
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .by= "site") 


discharge_all <- read_csv("data/stream_datasets/discharge_usgs.csv") %>% 
  mutate(date= as.Date(date_time)) %>% 
  summarise(discharge = mean(discharge), .by= c("site", "date")) %>% 
  left_join(sp_site_labels %>% select(siteID, USGSgageID), by = c("site" = "USGSgageID")) %>% 
  mutate(discharge_m3s = discharge*0.028316832)

sp_data <- read_csv("data/stream_datasets/stream_dataset.csv") %>% 
  mutate(date= as.Date(dateTimeUTC)) %>% 
  left_join(discharge_all, by = c("siteID", "date")) %>% 
  mutate(Discharge_joined = ifelse(is.na(Discharge_m3s) == TRUE, discharge_m3s, Discharge_m3s),
         siteID= str_remove(siteID, "/"))  %>% 
  left_join(sp_metab_data, by = c("siteID" = "site", "date"= "day_month"))
  


unique(sp_metab_data$site)

unique(discharge_all$siteID) 
unique(sp_data$siteID) 

#only ellipses, same axis
sp_data %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3))+
  stat_ellipse(geom = "polygon", alpha = .9, aes(fill = siteID))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  #geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  guides(fill = "none")+
  facet_wrap(~siteID)

ggsave("outputs/figures/streampulse_elipses.png", width = 12, height = 10)


sp_plot_points <- sp_data %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3))+
  stat_ellipse(geom = "polygon", alpha = .3, aes(fill = siteID))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  guides(fill = "none")+
  facet_wrap(~siteID, scales = "free")
  
ggsave("outputs/figures/streampulse_points.png", width = 12, height = 10)


sp_data %>% 
  drop_na(Discharge_joined) %>% 
  mutate(q.categories = case_when(Discharge_joined < 0.1 ~ "very small",
                                  Discharge_joined >= 0.1 & Discharge_joined < 1 ~ "small",
                                  Discharge_joined >= 1 & Discharge_joined < 10 ~ "medium",
                                  Discharge_joined >= 10 & Discharge_joined < 100 ~ "big",
                                  Discharge_joined >= 100 ~ "very big") %>% 
           as_factor() %>% 
           reorder(Discharge_joined, mean)) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3))+
  #geom_point(size=.1, alpha=.1, aes(color=q.categories))+
  stat_ellipse(geom = "polygon", alpha = .4, aes(colour=siteID, fill = q.categories), linewidth= 0)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  scale_fill_brewer(type= "seq", palette = "YlGnBu")+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)), title= "Each site one ellipse")+
  guides(colour = "none")

sp_data %>% 
  drop_na(Discharge_joined) %>% 
  mutate(q.categories = case_when(Discharge_joined < 0.1 ~ "very small",
                                  Discharge_joined >= 0.1 & Discharge_joined < 1 ~ "small",
                                  Discharge_joined >= 1 & Discharge_joined < 10 ~ "medium",
                                  Discharge_joined >= 10 & Discharge_joined < 100 ~ "big",
                                  Discharge_joined >= 100 ~ "very big") %>% 
           as_factor() %>% 
           reorder(Discharge_joined, mean)) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3))+
  #geom_point(size=.1, alpha=.1, aes(color=q.categories))+
  stat_ellipse(geom = "polygon", alpha = .4, aes(fill = q.categories), linewidth= 0)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  scale_fill_brewer(type= "seq", palette = "YlGnBu")+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)), title= "each ellipse one size bin")

sp_data %>% 
  filter(Discharge_joined> 0) %>% 
  summarise(across(c(Light_lux, GPP, WaterTemp_C, Discharge_joined, CO2dep_mmolm3, O2dep_mmolm3 ), 
                   list(mean= ~mean(., na.rm= TRUE), 
                        median= ~median(., na.rm= TRUE), 
                        cv= ~sd(., na.rm= TRUE)/mean(., na.rm= TRUE)*100)), .by = "siteID") %>% 
  ggplot(aes(GPP_mean, Discharge_joined_cv))+
  geom_point(aes(color=siteID, size= Discharge_joined_median))+
  geom_label(#data= . %>% filter(siteID %in% c("AbiskoM6", "SBM", "BRW", "BEC", "HBF")),
             aes(label= siteID), nudge_x=.5, nudge_y =.5)+
  labs(x= "Mean GPP (g O2 m-2 d-1)", y= "CV of discharge (%)")

#with ellipses by month
sp_data %>%
  mutate(month=month(dateTimeUTC), 
         season = case_when(month %in% c(12,1,2) ~ "winter",
                            month %in% 3:5 ~ "spring",
                            month %in% 6:8 ~ "summer",
                            month %in% 9:11 ~ "autumn")) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3, fill = month, group=month))+
  stat_ellipse(geom = "polygon", alpha = .3)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  #geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  scale_fill_viridis_c(option="inferno")+
  #guides(fill = "none")+
  facet_wrap(~siteID, scales = "free")

ggsave("outputs/figures/streampulse_elipses_month.png", width = 12, height = 10)

data <- sp_data %>% filter(siteID == "NHC") %>% 
  mutate(month = month(dateTimeUTC))
  


metrics_monthly <- sp_data %>% 
  drop_na(CO2dep_mmolm3, O2dep_mmolm3) %>% 
  mutate(year = year(dateTimeUTC),
         month = month(dateTimeUTC)) %>% 
  group_by(siteID, year, month) %>% 
  summarise(across(c(Discharge_joined, WaterTemp_C, Light_lux, Depth_m, pH, Nitrate_mgL, CDOM_ppb ), \(x) median(x, na.rm = TRUE)),
            meanCO2= mean(CO2dep_mmolm3),
            meanO2 = mean(O2dep_mmolm3), 
            n = n(),
            offset = meanCO2 + meanO2,
            elipse = elipse(CO2dep_mmolm3, O2dep_mmolm3),
            reg = regression_metrics(CO2dep_mmolm3, O2dep_mmolm3)  ) %>% 
  unnest(c(elipse, reg))

metrics_weekly <- sp_data %>% 
  drop_na(CO2dep_mmolm3, O2dep_mmolm3) %>% 
  mutate(year = year(dateTimeUTC),
         week = week(dateTimeUTC)) %>% 
  group_by(siteID, year, week) %>% 
  summarise(across(c(Discharge_joined, WaterTemp_C, Light_lux, Depth_m, pH, Nitrate_mgL, CDOM_ppb ), \(x) median(x, na.rm = TRUE)),
            meanCO2= mean(CO2dep_mmolm3),
            meanO2 = mean(O2dep_mmolm3), 
            n = n(),
            offset = meanCO2 + meanO2,
            elipse = elipse(CO2dep_mmolm3, O2dep_mmolm3),
            reg = regression_metrics(CO2dep_mmolm3, O2dep_mmolm3)  ) %>% 
  unnest(c(elipse, reg))



metrics_weekly %>% ungroup() %>% 
  summarise(across(everything(), \(x) median(x, na.rm = TRUE)), .by = "siteID") %>% 
ggplot(aes(Discharge_joined, meanCO2, color=WaterTemp_C ))+
 # geom_smooth()+
  geom_point()+ 
  scale_x_log10()


metrics_weekly %>% ungroup() %>% 
  summarise(across(everything(), \(x) median(x, na.rm = TRUE)), .by = "siteID") %>% 
  ggplot(aes(WaterTemp_C, meanCO2 ))+
  # geom_smooth()+
  geom_point()






selected_sites <- metrics_weekly %>% filter(siteID %in% c("AbiskoM6", "SBM", "BRW", "BEC"))


selected_sites %>% 
  mutate(site_type= case_when(siteID == "AbiskoM6"~ "AbiskoM6 - Dark & Stable",
                              siteID == "BEC"~ "Taconazo - Light & Satble",
                              siteID == "BRW"~ "BRW - Light & Stormy",
                              siteID == "SBM"~ "SBM - Dark & Stormy")) %>% 
  ggplot(aes(Discharge_joined,length, color= abs(corr)))+
  geom_point()+
  scale_color_viridis_c()+
  facet_wrap(~site_type, scales = "free_x")

sp_data %>%
  filter(siteID %in% c("AbiskoM6", "SBM", "BRW", "BEC")) %>% 
  mutate(month=month(dateTimeUTC), 
         season = case_when(month %in% c(12,1,2) ~ "winter",
                            month %in% 3:5 ~ "spring",
                            month %in% 6:8 ~ "summer",
                            month %in% 9:11 ~ "autumn"),
         site_type= case_when(siteID == "AbiskoM6"~ "AbiskoM6 - Dark & Stable",
                              siteID == "BEC"~ "Taconazo - Bright & Stable",
                              siteID == "BRW"~ "BRW - Bright & Stormy",
                              siteID == "SBM"~ "SBM - Dark & Stormy") %>% 
           as_factor() %>% 
           fct_relevel(c("AbiskoM6 - Dark & Stable", "Taconazo - Bright & Stable",
                         "SBM - Dark & Stormy", "BRW - Bright & Stormy" ))) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3, fill = month, group=month))+
  stat_ellipse(geom = "polygon", alpha = .3)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  #geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  scale_fill_viridis_c(option="inferno")+
  #guides(fill = "none")+
  facet_wrap(~site_type)



#old stuff



out <- sp_data %>% 
  drop_na(CO2dep_mmolm3, O2dep_mmolm3) %>% 
  mutate(year = year(dateTimeUTC),
         week = week(dateTimeUTC)) %>% 
  group_by(siteID, year, week) %>% 
  summarise(across(c(Discharge_m3s,WaterTemp_C, Light_lux, Depth_m, pH, Nitrate_mgL, CDOM_ppb ), \(x) median(x, na.rm = TRUE)),
            meanCO2= mean(CO2dep_mmolm3),
            meanO2 = mean(O2dep_mmolm3), 
            n = n(),
            offset = meanCO2 + meanO2,
            elipse = elipse(CO2dep_mmolm3, O2dep_mmolm3),
            reg = regression_metrics(CO2dep_mmolm3, O2dep_mmolm3)  ) %>% 
  unnest(c(elipse, reg))
  
unique(out$siteID)

out %>% 
  ungroup(siteID, year, week) %>% 
  mutate(date_week = as.Date(paste(year, week-1, 1), format = "%Y %U %u")) %>% 
  #filter(siteID == "BEC") %>% 
  filter(Discharge_m3s> 0, EQ < 2) %>% 
  pivot_longer(offset:EQ, names_to = "metric", values_to = "value") %>% 
  ggplot(aes(Discharge_m3s, value, color= abs(corr)))+
  geom_point()+
  scale_color_viridis_c()+
  facet_wrap(~metric, scales = "free")

out %>% 
  ungroup(siteID, year, week) %>% 
  mutate(date_week = as.Date(paste(year, week-1, 1), format = "%Y %U %u")) %>% 
  filter(siteID == "BEC") %>% 
  filter(EQ < 2) %>% 
  pivot_longer(offset:EQ, names_to = "metric", values_to = "value") %>% 
  ggplot(aes(WaterTemp_C, value, color= abs(corr)))+
  geom_point()+
  scale_color_viridis_c()+
  facet_wrap(~metric, scales = "free")

out %>% 
  ungroup(siteID, year, week) %>% 
  filter(Discharge_m3s> 0) %>% 
  mutate(date_week = as.Date(paste(year, week-1, 1), format = "%Y %U %u")) %>% 
  ggplot(aes(Discharge_m3s, offset, color=WaterTemp_C ))+
  geom_point()+
  scale_color_viridis_c()+
  facet_wrap(~siteID, scales = "free")


out %>% 
  ungroup(siteID, year, week) %>% 
 #filter(abs(corr) > .7) %>% 
  mutate(date_week = as.Date(paste(year, week-1, 1), format = "%Y %U %u")) %>% 
  ggplot(aes(date_week, length))+
  geom_point()+
  #scale_color_viridis_c()+
  facet_wrap(~siteID, scales = "free")

out %>% 
  ungroup(siteID, year, week) %>% 
  select(siteID, offset:EQ) %>% 
  pivot_longer(width:EQ, names_to = "metric", values_to = "value") %>% 
  ggplot(aes(offset, value, color= siteID))+
  geom_point()+
  facet_wrap(~metric, scales = "free")


out %>% 
  ungroup(siteID, year, week) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .by = siteID) %>% 
  left_join(sp_site_data) %>% 
  select(where(is.numeric)) %>%  
  correlate() %>% 
  focus(length) %>% 
  arrange(desc(abs(length))) %>% 
  print(n = 50)


out %>% 
  ungroup(siteID, year, week) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .by = siteID) %>% 
  left_join(sp_site_data) %>% 
  ggplot(aes(Discharge_m3s, offset))+
  geom_point()+
  scale_x_log10()


out %>% 
  ungroup(siteID, year, week) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .by = siteID) %>% 
  left_join(sp_site_data) %>% 
  ggplot(aes(crop_cover_sub_per, offset))+
  geom_point()

out %>% 
  ungroup(siteID, year, week) %>% 
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .by = siteID) %>% 
  left_join(sp_site_data) %>% 
  ggplot(aes(crop_cover_up_per, length))+
  geom_point()


#with ellipses by month
sp_data %>%
  mutate(month=month(dateTimeUTC), 
         season = case_when(month %in% c(12,1,2) ~ "winter",
                            month %in% 3:5 ~ "spring",
                            month %in% 6:8 ~ "summer",
                            month %in% 9:11 ~ "autumn")) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3, fill = month, group=month))+
  stat_ellipse(geom = "polygon", alpha = .3)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  #geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  scale_fill_viridis_c(option="inferno")+
  #guides(fill = "none")+
  facet_wrap(~siteID, scales = "free")


#with ellipses by month
example_data <- sp_data %>%
  filter(siteID %in% c("BRW", "HBF", "Abisko/M10", "Drain")) %>% 
  mutate(month=month(dateTimeUTC), 
         season = case_when(month %in% c(12,1,2) ~ "winter",
                            month %in% 3:5 ~ "spring",
                            month %in% 6:8 ~ "summer",
                            month %in% 9:11 ~ "autumn"),
         site2 = case_when(siteID == "BRW" ~ "site1",
                           siteID == "HBF" ~ "site2",
                           siteID == "Abisko/M10" ~ "site3",
                           siteID == "Drain" ~ "site4"))  
example_data %>% 
ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3, fill = month, group=month))+
 # geom_point(aes(color=season), size=.1, alpha=.1)+
  stat_ellipse(geom = "polygon", alpha = .3)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  #geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  scale_fill_viridis_c(option="inferno")+
  facet_wrap(~site2, scales = "free")

example_data %>%  
  filter(site2 == "site4") %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3))+
  # geom_point(aes(color=season), size=.1, alpha=.1)+
  stat_ellipse(geom = "polygon", alpha = .3)+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  #geom_point(size=.1, alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mmol~m^-3)), y=expression(O[2]~departure~(mmol~m^-3)))+
  #scale_fill_viridis_c(option="inferno")+
  facet_wrap(~month)+
  labs(title= "site 4")

data <- sp_data %>% filter(siteID == "NHC") %>% 
  mutate(month = month(dateTimeUTC))



out <- sp_data %>% 
  drop_na(CO2dep_mmolm3, O2dep_mmolm3) %>% 
  mutate(year = year(dateTimeUTC),
         week = week(dateTimeUTC)) %>% 
  group_by(siteID) %>% 
  summarise(across(c(Discharge_m3s,WaterTemp_C, Light_lux, Depth_m, pH, Nitrate_mgL, CDOM_ppb ), \(x) mean(x, na.rm = TRUE)),
            meanCO2= mean(CO2dep_mmolm3),
            meanO2 = mean(O2dep_mmolm3), 
            n = n(),
            offset = meanCO2 + meanO2,
            elipse = elipse(CO2dep_mmolm3, O2dep_mmolm3),
            reg = regression_metrics(CO2dep_mmolm3, O2dep_mmolm3)  ) %>% 
  unnest(c(elipse, reg))

out %>% left_join(sp_site_data) %>% 
  ggplot(aes(CDOM_ppb, offset))+
  geom_point(size= 3)

library(corrr)

out %>% 
  left_join(sp_site_data) %>% 
  select(where(is.numeric)) %>% 
  filter(Discharge_m3s > 0) %>% 
  correlate() %>% 
  focus(offset) %>% 
  arrange(desc(abs(offset))) %>% 
  print(n=40)

out %>% 
  left_join(sp_site_data) %>% 
  select(where(is.numeric)) %>% 
  filter(Discharge_m3s > 0) %>% 
  correlate() %>% 
  focus(width) %>% 
  arrange(desc(abs(width))) %>% 
  print(n=20)


out %>% 
  left_join(sp_site_data) %>% 
  ggplot(aes(wetland_up_per_u01, offset))+
  geom_point()


out %>% 
  left_join(sp_site_data) %>% 
  ggplot(aes(slope_up_degrees, width))+
  geom_point(size=3)


