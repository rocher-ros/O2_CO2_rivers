# Load  and install packages ----
# List of all packages needed
package_list <- c('tidyverse', 'ellipse', 'lmodel2', 'e1071', 'patchwork')

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
sp_site_labels <- read_csv("prepared data/river data/streampulse/streampulse_site_data.csv") %>% 
  filter(str_detect(variableList,"CO2_ppm" ))

sp_site_data <- read_csv("prepared data/river data/streampulse/site_info_streams.csv") %>% 
  select(-geometry)


sp_metab_data <- read_csv("prepared data/river data/streampulse/all_daily_model_results.csv") %>% 
  filter(GPP>= 0, ER <= 0, valid_day == 1) %>% 
  mutate(day_month= as.Date(date))

sp_metab_avg <- sp_metab_data %>% 
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .by= "site") 


discharge_all <- read_csv("prepared data/river data/streampulse/discharge_usgs.csv") %>% 
  mutate(date= as.Date(date_time)) %>% 
  summarise(discharge = mean(discharge), .by= c("site", "date")) %>% 
  left_join(sp_site_labels %>% select(siteID, USGSgageID), by = c("site" = "USGSgageID")) %>% 
  mutate(discharge_m3s = discharge*0.028316832)

sp_data <- read_csv("prepared data/river data/streampulse/stream_dataset.csv") %>% 
  mutate(date= as.Date(dateTimeUTC)) %>% 
  left_join(discharge_all, by = c("siteID", "date")) %>% 
  mutate(Discharge_joined = ifelse(is.na(Discharge_m3s) == TRUE, discharge_m3s, Discharge_m3s),
         siteID= str_remove(siteID, "/"))  %>% 
  left_join(sp_metab_data, by = c("siteID" = "site", "date"= "day_month")) %>% 
  mutate(siteID = case_when(siteID == "k4" ~ "KC4",
                            siteID == "k6" ~ "KC6",
                            siteID == "k7" ~ "KC7",
                            TRUE ~siteID),
    GPP= ifelse(siteID %in%  c("KC4", "KC6", "KC7"), 0.02, GPP),
    ER= ifelse(siteID %in%  "KC6", -4, ER))
  


metab_q_avg_metrics <- sp_data %>% 
  filter(Discharge_joined> 0) %>% 
  summarise(across(c(GPP, ER, Discharge_joined ), 
                   list(mean= ~mean(., na.rm= TRUE), 
                        median= ~median(., na.rm= TRUE), 
                        cv= ~sd(., na.rm= TRUE)/mean(., na.rm= TRUE)*100,
                        skewness= ~skewness(., na.rm= TRUE)
                        )), .by = "siteID")  %>% 
  mutate(discharge_category= case_when(Discharge_joined_median <= 0.1 ~"[0 - 0.1]",
                                       Discharge_joined_median > 0.1 & Discharge_joined_median <= 1  ~"(0.1 - 0.5]",
                                       Discharge_joined_median > 1 & Discharge_joined_median <= 5  ~"(1 - 5]",
                                       Discharge_joined_median > 5 & Discharge_joined_median <= 10  ~"(5 - 10]",
                                       Discharge_joined_median > 10 ~"(10 - 100]") %>% 
           as_factor() %>% fct_reorder(Discharge_joined_median))
  

metab_q_avg_metrics %>% 
  ggplot(aes(GPP_mean, Discharge_joined_skewness))+
  geom_point(aes(color=discharge_category), size= 3 )+
  geom_label(data= . %>% filter(siteID %in% c("KC6","SBM", "BRW", "BEC")),
             aes(label= siteID), nudge_x=.3, label.size  = NA, alpha = 0)+
  labs(x= expression(GPP~(g~O[2]~m^-2~d^-1)), y= "Skewness of discharge", color= expression((Discharge~m^3~s^-1)))+
  scale_color_viridis_d(direction = -1, option= "mako", begin= .1, end= .9)+
  theme_classic()+
  theme(legend.position = c(.85,.75))


ggsave("plots/SM/fig_gpp_sqew.png", width = 6, height = 4)

metab_q_avg_metrics %>% filter(siteID %in% c("KC6","SBM", "BRW", "BEC")) %>% 
  select(siteID, GPP_mean, ER_mean, Discharge_joined_mean, Discharge_joined_skewness)

                               
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




selected_sites <-  c("KC6", "SBM", "BRW", "BEC")

metrics_site <- sp_data %>% 
  filter(siteID %in% selected_sites) %>% 
  drop_na(CO2dep_mmolm3, O2dep_mmolm3) %>% 
  group_by(siteID) %>% 
  summarise(meanCO2= mean(CO2dep_mmolm3),
            meanO2 = mean(O2dep_mmolm3), 
            n = n(),
            offset = meanCO2 + meanO2,
            elipse = elipse(CO2dep_mmolm3, O2dep_mmolm3),
            reg = regression_metrics(CO2dep_mmolm3, O2dep_mmolm3)  ) %>% 
  unnest(c(elipse, reg)) %>% 
  mutate(site_type= case_when(siteID == "KC6"~ "Dark and Stable",
                              siteID == "BEC"~ "Bright and Stable",
                              siteID == "BRW"~ "Bright and Stormy",
                              siteID == "SBM"~ "Dark and Stormy") %>% 
           as_factor() %>% 
           fct_relevel(c("Dark and Stable", "Bright and Stable",
                         "Dark and Stormy", "Bright and Stormy" )))


metrics_plot <- 
  metrics_site %>% 
  mutate(text.string = paste0("Centroid: (",round(meanCO2,0), ";", round(meanO2,0), ")","\n",
                              "Length: ", round(length, 0),"\n",
                              "Offset: ", round(offset, 0),"\n",
                              "EQ: ", round(EQ,1), "\n",
                              expression(R^2),": ", round(corr^2,2)),
         text.site = c("Black Earth creek (WI, USA)", "Brewery creek (WI, USA)", "Krycklan C6 (Sweden)", "Saddleback mountain (NH, USA)") )

metrics_plot %>% select(siteID, text.string)


sp_data %>%
  filter(siteID %in% selected_sites) %>% 
  mutate(site_type= case_when(siteID == "KC6"~ "Dark and Stable",
                              siteID == "BEC"~ "Bright and Stable",
                              siteID == "BRW"~ "Bright and Stormy",
                              siteID == "SBM"~ "Dark and Stormy") %>% 
           as_factor() %>% 
           fct_relevel(c("Dark and Stable", "Bright and Stable",
                         "Dark and Stormy", "Bright and Stormy" ))) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3,  color= site_type, group=site_type))+
  #geom_point(size=.05, alpha=.05)+
  stat_ellipse( geom = "polygon", alpha = 0, level = 0.9)+
  geom_text(data= metrics_plot, aes(x=610, y= 190, label=text.site), size= 3, color= "black", hjust = 1)+
  stat_density_2d_filled(aes(alpha = after_stat(level), fill= site_type), linewidth= .1,
                         geom = "polygon", contour = TRUE, contour_var = "ndensity",
                         breaks = c(0.02, 0.05, 0.25, 0.5,0.75,1))+
  scale_fill_manual(values = c("cadetblue4", "darkolivegreen4", "coral4", "goldenrod3"))+
  scale_color_manual(values = c("cadetblue4", "darkolivegreen4", "coral4", "goldenrod3"))+
  geom_hline(yintercept = 0, linetype=1)+
  geom_vline(xintercept = 0, linetype=1)+
  geom_abline(slope=-1, intercept = 0, linetype=2)+
  labs(x=expression(CO[2]~departure~(mu*mol~L^-1)), y=expression(O[2]~departure~(mu*mol~L^-1)))+
  facet_wrap(~site_type)+
  theme_classic()+
  theme(legend.position ="none",
        panel.border =element_rect(linewidth = 1, fill= NA),
        strip.text = element_text(size=12, face= "bold"),
        strip.background = element_blank()) 

ggsave("plots/main/fig3_met_regimes.png", width = 7, height = 6)


