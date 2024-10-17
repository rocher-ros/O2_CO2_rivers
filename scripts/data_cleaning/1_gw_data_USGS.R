
# Load packages
#library(dataRetrieval)
library(tidyverse)
library(janitor)

resultphyschem <- read_csv("empirical data/groundwater data/_raw/resultphyschem.csv")

sites <- read_csv("empirical data/groundwater data/_raw/station.csv")


names(sites)
names(resultphyschem)

unique(resultphyschem$CharacteristicName) 

unique(sites$MonitoringLocationTypeName)

chem_sites <- resultphyschem %>% 
  left_join(sites, by= "MonitoringLocationIdentifier") %>% 
  filter(MonitoringLocationTypeName %in% c("Subsurface", "Subsurface: Unsaturated zone",
                                           "Subsurface: Groundwater drain", "Aggregate groundwater use"))


chem_site_avg <- chem_sites %>% 
  drop_na(`ResultMeasure/MeasureUnitCode`) %>% 
  mutate(variable_unit = paste0(CharacteristicName,"_",`ResultMeasure/MeasureUnitCode`)) %>% 
  pivot_wider(names_from =variable_unit, values_from =ResultMeasureValue ) %>% 
  summarise(across(`pH_std units`:`Inorganic carbon_%`, ~mean(.x, na.rm = TRUE)), 
            latitude = mean(LatitudeMeasure, na.rm = TRUE),
            longitude = mean(LongitudeMeasure,  na.rm = TRUE),
            depth = mean(`VerticalMeasure/MeasureValue`, na.rm = TRUE),
            across(where(is.character), first),
            .by= "MonitoringLocationIdentifier") %>% 
  clean_names()  %>% 
  mutate(carbon_dioxide_mg_l = ifelse(carbon_dioxide_mg_l > 150, NA, carbon_dioxide_mg_l),#extreme pco2 values, romving above 150.000 ppm
         co2_umol_l = carbon_dioxide_mg_l/44*1e+3,
         co2_ppm = co2_umol_l*1e+3* 0.04477565,
         oxygen_mg_l = ifelse(oxygen_mg_l > 13, NA, oxygen_mg_l))

names(chem_site_avg)


write_csv(chem_site_avg, "empirical data/groundwater data/usgs_gw_co2_o2_site_avg.csv")

chem_site_avg %>% drop_na(oxygen_mg_l) %>% 
ggplot( aes(longitude, latitude, color= oxygen_mg_l))+
  geom_point()

chem_site_avg %>% filter(co2_ppm < 180000) %>% 
  ggplot( aes(longitude, latitude, color= co2_ppm))+
  geom_point()

chem_site_avg %>% 
  ggplot( aes(carbon_dioxide_mg_l, co2_ppm))+
  geom_point()






