

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)


#groundwater data 

gw_usgs <- read_csv("empirical data/groundwater data/usgs_gw_co2_o2_site_avg.csv") %>% 
  mutate(carbon_dioxide_mg_l = ifelse(carbon_dioxide_mg_l > 150, NA, carbon_dioxide_mg_l),#extreme pco2 values, romving above 150.000 ppm
         co2_umol_l = carbon_dioxide_mg_l/44*1e+3,
         co2_ppm = co2_umol_l*1e+3* 0.04477565,
         oxygen_mg_l = ifelse(oxygen_mg_l > 13, NA, oxygen_mg_l),
         oxygen_umol_l= oxygen_mg_l/32*1e+3) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

gw_plotting <- gw_usgs %>% mutate(shape_plot= case_when(is.na(oxygen_umol_l) == TRUE & is.na(co2_umol_l) == TRUE ~ "none" ,
                                                        is.na(oxygen_umol_l) == FALSE & is.na(co2_umol_l) == TRUE ~ "Only CO2",
                                                        is.na(oxygen_umol_l) == TRUE & is.na(co2_umol_l) == FALSE ~ "Only O2",
                                                        is.na(oxygen_umol_l) == FALSE & is.na(co2_umol_l) == FALSE ~ "Both",
                                                        ) )
#sites form Appling 2019
daily_metab <- read_delim("empirical data/river data/Appling2019/daily_predictions.tsv") %>% 
  filter(GPP>=0, ER <= 0, GPP.Rhat < 1.1, ER.Rhat < 1.1, K600.Rhat < 1.1)

#get site average data for all variables, calculated in script # 4
site_metab <- read_csv("empirical data/river data/Appling2019/site_avgs_gwinputs.tsv") %>% 
  select(site_name, ends_with("_median"), gw_frac, spQ_mmday, tot_area) %>% 
  rename_with(~ str_remove(., "_median"), everything()) %>% 
  mutate(gw_frac = ifelse(gw_frac > 100 | gw_frac < 0, NA, gw_frac))

site_catchments <- read_sf("empirical data/river data/Appling2019/catchment_shapefile/catchment_shapefile.shp")
outlet_catchments <- read_sf("empirical data/river data/Appling2019/points_shapefile/points_shapefile.shp")

#read the water chemistry file, to get alkalinity
chemistry_usgs <- read_csv("empirical data/river data/USGS_data/chemistry_site_avg.csv")

#join to the previous file
met_with_chem <- site_metab %>% 
  mutate(site = str_remove_all(site_name, "nwis_") %>% 
           paste0("USGS-", .)) %>% 
  left_join(chemistry_usgs, by = "site") %>% 
  mutate(alk_mmol_m3 = alkalinity_mg_l_ca_co3*0.02*1000)


north_america <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(continent== "North America")

rivers50 <- ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")

plot_gw <- 
  ggplot()+
  geom_sf(data= north_america, linewidth= 0.1)+
  geom_sf(data= rivers50, color= "cornflowerblue")+
  #geom_sf(data= gw_usgs %>% drop_na(oxygen_umol_l), shape= 2)+
  #geom_sf(data= gw_usgs %>% drop_na(co2_umol_l), shape= 6)+
  geom_sf(data = gw_plotting %>% filter(shape_plot != "none"), aes(shape= shape_plot))+
  scale_shape_manual(values= c(11,6,2))+
  coord_sf(xlim = c(-170,-52), ylim = c(15,70), expand = FALSE)+
  theme_minimal()+
  theme(legend.position = c(0.92, 0.82), legend.box.background = element_rect(fill= "white"))+
  labs(title= "Groundwater sites", shape= "Groundwater")


plot_rivers <- ggplot()+
  geom_sf(data= north_america, linewidth= 0.1)+
  geom_sf(data= rivers50, color= "cornflowerblue")+
  geom_sf(data= outlet_catchments)+
  coord_sf(xlim = c(-170,-52), ylim = c(15,70), expand = FALSE)+
  theme_minimal()+
  labs(title= "River sites")


plot_rivers + plot_gw + plot_layout(ncol = 1) + plot_annotation(tag_levels = "a")

ggsave("plots/SM/maps_sites.png")


## Figure of alkalinity with river discharge 

met_with_chem %>% 
  ggplot(aes(discharge, alk_mmol_m3 ))+
  geom_smooth(color= "gray20", linetype= 2, alpha=.2)+
  geom_point()+
  scale_x_log10(labels = scales::number, breaks= c(0.01,0.1, 1, 10, 100, 1000))+
  theme_classic()+
  labs(x=expression(Discharge~(m^3~s^-1)), y= expression(Alkalinity~(mu*mol~L^-1)))

ggsave("plots/SM/alkalinity_discharge.png", width = 7, height = 5)
