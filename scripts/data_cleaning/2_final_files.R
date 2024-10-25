# Information ----
# Script to pre-process river GIS data and the USGS river chemistry data for later use.
# Author: Gerard Rocher-Ros
# Last edit: 2024-09-03

# Install and Load libraries ----

# List of all packages needed
package_list <- c('tidyverse', 'sf', 'nhdplusTools', 'janitor')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

# Now load all the packages
lapply(package_list, require, character.only = TRUE)


# 1. Download stream network properties from nhdplus ----

outlet_catchments <- read_sf("processed data/river data/Appling2019/points_shapefile/points_shapefile.shp")

data_nhd <- tibble(site_nm= outlet_catchments$site_nm,
                   tot_length= NA,
                   tot_area= NA,
                   slope= NA,
                   order = NA,
                   avg_dicharge = NA )


for(i in seq_along(data_nhd$site_nm)){
  
  start_point <- outlet_catchments[i,]
  
  start_comid <- discover_nhdplus_id(start_point)
  
  flowline <- navigate_nldi(list(featureSource = "comid", 
                                 featureID = start_comid), 
                            mode = "upstreamTributaries", 
                            distance_km = 10)
  
  subset_file <- tempfile(fileext = ".gpkg")
  subset <- subset_nhdplus(comids = as.integer(flowline$UT$nhdplus_comid),
                           output_file = subset_file,
                           nhdplus_data = "download", 
                           flowline_only = TRUE,
                           return_data = TRUE, overwrite = TRUE)
  
  flowline <- subset$NHDFlowline_Network
  
  names(flowline)  
  data_nhd$tot_length[i] = flowline$arbolatesu[flowline$comid == start_comid]
  data_nhd$tot_area[i] = flowline$totdasqkm[flowline$comid == start_comid]
  data_nhd$slope[i]= flowline$slope[flowline$comid == start_comid]
  data_nhd$order[i]= flowline$streamorde[flowline$comid == start_comid]
  data_nhd$avg_dicharge[i]= flowline$qa_ma[flowline$comid == start_comid]
  
  print(paste("done", i, "of 602"))
}


#export the file for further analysis
write_csv(data_nhd, "processed data/river data/USGS_data/catchment_lengths_properties.csv")


# 2. Process groundwater files ----
#read files downloaded from USGS waterdata portal
resultphyschem <- read_csv("raw data/groundwater data/_raw/resultphyschem.csv")

sites <- read_csv("raw data/groundwater data/_raw/station.csv")

# check file names
names(sites)
names(resultphyschem)

unique(resultphyschem$CharacteristicName) 

unique(sites$MonitoringLocationTypeName)

#remove sites from mines and caves
chem_sites <- resultphyschem %>% 
  left_join(sites, by= "MonitoringLocationIdentifier") %>% 
  filter(MonitoringLocationTypeName %in% c("Subsurface", "Subsurface: Unsaturated zone",
                                           "Subsurface: Groundwater drain", "Aggregate groundwater use"))

#wrangle the file, outliers and so forth
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
  mutate(carbon_dioxide_mg_l = ifelse(carbon_dioxide_mg_l > 150, NA, carbon_dioxide_mg_l),#extreme pco2 values, romving above 150 mgL
         co2_umol_l = carbon_dioxide_mg_l/44*1e+3,
         co2_ppm = co2_umol_l*1e+3* 0.04477565,
         oxygen_mg_l = ifelse(oxygen_mg_l > 13, NA, oxygen_mg_l))

names(chem_site_avg)


write_csv(chem_site_avg, "processed data/groundwater data/usgs_gw_co2_o2_site_avg.csv")


# 3. process water chemistry data for all sites 
chemistry_all <- list.files(path = "raw data/river data/USGS_data/nwis/",
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

vars_all <- chemistry_all %>% summarise(n= n(), .by = "variable")

vars_all %>% arrange(desc(n)) %>% #filter(str_detect(variable,"Chloro") == TRUE) %>% 
  print(n=1000) 

#select some potentially useful variables for the study
useful_vars <- c("Chemical oxygen demand, (low level)" , "Nitrate", "Chloride", "Chlorophyll a", "Inorganic nitrogen (nitrate and nitrite)",
                 "Phosphorus", "Strontium", "Carbonate", "Kjeldahl nitrogen", "Iron" , "Calcium" , "Sodium", "Oxygen", "Inorganic carbon",
                 "Bicarbonate", "Sulfate", "Acidity, (H+)", "pH", "Stream flow, instantaneous" , "Temperature, water",
                 "Color", "Alkalinity", "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)", "Suspended Sediment Concentration (SSC)",
                 "Turbidity", "Biomass, periphyton", "Phytoplankton", "Specific conductance", "Carbon dioxide", "Orthophosphate",
                 "Organic carbon", "Temperature, air, deg C", "Stream flow, mean. daily" , "Barometric pressure", "Depth",
                 "Organic Nitrogen" , "Chemical oxygen demand", "Biochemical oxygen demand, standard conditions", "UV 254", "Caffeine")

#have a quick look
chemistry_all %>% 
  filter(variable == "Chlorophyll a")

#keep only water  samples
chemistry_clean <- chemistry_all %>% 
  filter(variable %in% useful_vars, variable.info != c("Bed Sediment", "Recoverable")) %>% 
  mutate(value = as.numeric(value)) 

#do site averages and export it
avg_site <- chemistry_clean %>% 
  summarise(value_mean = mean(value , na.rm = TRUE),
            unit= first(unit),
            .by= c("site", "variable")) %>% 
  pivot_wider(names_from = c(variable, unit),
              values_from= value_mean) %>% 
  janitor::clean_names()

write_csv(avg_site, "processed data/river data/USGS_data/chemistry_site_avg.csv")







