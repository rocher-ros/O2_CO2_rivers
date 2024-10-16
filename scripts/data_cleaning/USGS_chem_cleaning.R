
library(tidyverse)

chemistry_all <- list.files(path = "model/data/nwis/raw_data_nwis/",
                            pattern="*.csv", 
                            full.names = T) %>% 
  map_df(~read_csv(., col_types = cols(.default = col_character()))) %>%  
  mutate(date= as.character(ActivityStartDate),
         time= as.character(ActivityStartTime.Time),
         site= as.character(MonitoringLocationIdentifier),
         conditions = as.character(HydrologicCondition),
         variable = as.character(CharacteristicName),
         variable.info= as.character(ResultSampleFractionText),
         value = as.character(ResultMeasureValue),
         unit = as.character(ResultMeasure.MeasureUnitCode), 
         .keep = "none")

options(max.print=2000)
unique(chemistry_all$variable) 

vars_all <- chemistry_all %>% summarise(n= n(), .by= "variable")

vars_all %>% arrange(desc(n)) %>% #filter(str_detect(variable,"Chloro") == TRUE) %>% 
  print(n=1000) 

useful_vars <- c("Chemical oxygen demand, (low level)" , "Nitrate", "Chloride", "Chlorophyll a", "Inorganic nitrogen (nitrate and nitrite)",
                 "Phosphorus", "Strontium", "Carbonate", "Kjeldahl nitrogen", "Iron" , "Calcium" , "Sodium", "Oxygen", "Inorganic carbon",
                 "Bicarbonate", "Sulfate", "Acidity, (H+)", "pH", "Stream flow, instantaneous" , "Temperature, water",
                 "Color", "Alkalinity", "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)", "Suspended Sediment Concentration (SSC)",
                 "Turbidity", "Biomass, periphyton", "Phytoplankton", "Specific conductance", "Carbon dioxide", "Orthophosphate",
                 "Organic carbon", "Temperature, air, deg C", "Stream flow, mean. daily" , "Barometric pressure", "Depth",
                 "Organic Nitrogen" , "Chemical oxygen demand", "Biochemical oxygen demand, standard conditions", "UV 254", "Caffeine")

chemistry_all %>% 
  filter(variable == "Chlorophyll a")

chemistry_clean <- chemistry_all %>% 
  filter(variable %in% useful_vars, variable.info != c("Bed Sediment", "Recoverable")) %>% 
  mutate(value = as.numeric(value)) 

avg_site <- chemistry_clean %>% 
  summarise(value_mean = mean(value , na.rm = TRUE),
            unit= first(unit),
            .by= c("site", "variable")) %>% 
  pivot_wider(names_from = c(variable, unit),
              values_from= value_mean) %>% 
  janitor::clean_names()

write_csv(avg_site, "empirical data/river data/USGS_data/chemistry_site_avg.csv")
