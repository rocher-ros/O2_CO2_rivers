# Custom functions for CO2 and O2 ----

#Calculate CO2 from ppm to mmolm3
ppm_to_mmolm3 <- function(CO2_ppm,atm_pressure,temp.water ) { #atm_pressure in atmospheres
  Hconstant <- exp(-58.0931+(90.5069*(100/(temp.water+273.15)))+(22.294*log((temp.water+273.15)/100))) # [mol dm-3 atm-1]
  
  CO2_ppm*atm_pressure*Hconstant#in mmol CO2 /m3
}

# OSAT FUNCTION 
#oxygen saturation function based on temp, salinity and pressure as in Grace et al., 2015 L&O
#pressure in atmospheres and temp in celsius
osat_grace <- function(temp, bp, sal){ 
  t_K = 273.15 + temp
  S1 = 157570.1/t_K
  S2 = -66423080/(t_K^2)
  S3 = 12438000000/(t_K^3)
  S4 = -862194900000/(t_K^4)
  S5 = -1*sal*(0.017674 - (10.754/t_K) + (2140.7/(t_K^2)))
  
  f_sal = exp(-139.34411 + S1 +S2 + S3 + S4 + S5)
  
  alpha = 0.00975 - 0.0001426*t_K + 0.0000006436*t_K^2
  beta = exp( 11.8571 - 3840.77/t_K - 216961/(t_K^2) )
  gamma_p = ((1-beta/bp)/(1-beta) ) * ( (1-alpha*bp)/(1-alpha) )
  sato = f_sal*bp*gamma_p
  
  sato
}

#estimate atompsheric pressure (in atm) from elevation (m)
estimate_pressure <- function(elevation){
  pressure <- (1-(.0000225577*elevation))^5.25588
  return(pressure)
}

#read many files and keep name as column
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(siteID = flnm)
}

# Load  and install packages ----
# List of all packages needed
package_list <- c('tidyverse', "elevatr", "lubridate", "sf", "readxl")

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

# Now load all the packages
lapply(package_list, require, character.only = TRUE)


# Define constants
air_CO2 = 400 #ppm

# Process streampulse data ----
# Read files 
sp_data <- read_csv("empirical data/river data/streampulse/_raw/streampulse/streampulse_sites_with_co2_wide.csv") %>% 
  drop_na(CO2_ppm, DO_mgL, WaterTemp_C) %>%
  filter(CO2_ppm > 100, CO2_ppm < 20000, WaterTemp_C > -10, WaterTemp_C < 40, DO_mgL > 1,
         !siteID %in% c("AbiskoM1", "AbiskoM6","AbiskoM9", "AbiskoM10","AbiskoM17", "AbiskoM16")) %>% 
  mutate(AirPres_atm = AirPres_kPa/101.325)

#read the new ones form Miellajokka to replace in the streampulse DB
miellajokka <- list.files("empirical data/river data/streampulse/_raw/miellajokka/", full.names = TRUE) %>% 
  map_df(~read_plus(.)) %>% 
  mutate(siteID = str_replace(siteID,"empirical data/river data/streampulse/_raw/miellajokka/", "Abisko"),
         siteID = str_remove(siteID,"_2016_co2.csv"),
         siteID = str_remove(siteID,"_2015_co2.csv"),
         siteID = str_remove(siteID,"/"),
         dateTimeUTC = local.time,
         SpecCond_uScm = NA,
         WaterTemp_C = temp.water,
         Depth_m = depth1,
         Light_lux = NA,
         Discharge_m3s = discharge.daily/1000,
         AirPressCombined_atm = airpress/101.325,
         CO2_mmolm3 = ppm_to_mmolm3(co2_ppm, AirPressCombined_atm, WaterTemp_C),
         O2_mmolm3 = DO.obs/32*1000,
         CO2eq_mmolm3 = ppm_to_mmolm3(air_CO2, AirPressCombined_atm, WaterTemp_C),
         O2eq_mmolm3 = osat_grace(WaterTemp_C, AirPressCombined_atm, sal =0.001)/32*1000,
         CO2dep_mmolm3 = CO2_mmolm3 - CO2eq_mmolm3,
         O2dep_mmolm3 = O2_mmolm3 - O2eq_mmolm3,
         CDOM_ppb=NA,
         pH=NA, 
         Turbidity_NTU=NA, 
         Turbidity_FNU=NA, 
         SpecCond_uScm=NA, 
         Nitrate_mgL=NA ) %>% 
  select(siteID, dateTimeUTC, WaterTemp_C, Depth_m,  Light_lux, Discharge_m3s, AirPressCombined_atm, CO2_mmolm3:O2dep_mmolm3, 
         CDOM_ppb, pH, Turbidity_NTU, Turbidity_FNU, SpecCond_uScm, Nitrate_mgL)


sp_site_data <- read_csv("empirical data/river data/streampulse/_raw/streampulse/streampulse_site_data.csv") %>% 
  filter(str_detect(variableList,"CO2_ppm" )) 


sites_sf <- sp_site_data %>% 
  st_as_sf( coords = c("longitude", "latitude"),  crs = 4326) 


for(i in seq_along(sites_sf$siteID)){
  sites_sf$elev[i] <- get_elev_point(locations = sites_sf[i,], prj = st_crs(sites_sf), src = "aws", z = 10)$elevation 
  }
  


#Add elevation manually for some sites 
sp_site_data <- sp_site_data %>% 
  left_join(sites_sf %>% st_drop_geometry() %>% select(siteID, elev)) %>% 
  mutate( AirPresElev_atm = estimate_pressure(elev) )


sp_site_data %>% 
  filter(str_detect(variableList,"CO2_ppm" )) %>% 
  select(siteID) %>% 
  print(n=35)


#get the kyrcklan DO data form Gomez-Gener et al. 2021 NatComm
krycklan_do <- read_csv("empirical data/river data/streampulse/_raw/krycklan/raw/gomez_gener_ts_DO_Krycklan.csv")




#read in the co2 data form kryclan, obtained from monitoring
c7_co2 <- read_csv("empirical data/river data/streampulse/_raw/krycklan/raw/C7__Q_TW_CO2_2017-2022.Time_Series_Data.2024030405212542.csv", 
                skip = 13) %>% 
  rename(co2_ppm = `CO2 (Dis)@110107_C7`, q_ls = `Discharge@110107_C7`,
        temp_water_c = `Water Temp@110107_C7`, datetime = TimeStamp) %>% 
  mutate(year= year(datetime), month= month(datetime), day= day(datetime), hour= hour(datetime),
         date_hour = make_datetime(year, month, day, hour, 0, 0)) %>% 
  select(date_hour, co2_ppm, q_ls, temp_water_c) %>% 
  summarise(across(everything(), mean), .by ="date_hour") %>% 
  mutate(site= "k7")

c18_co2 <- read_csv("empirical data/river data/streampulse/_raw/krycklan/raw/C18_Q_TW_CO2_2017-2022.Time_Series_Data.2024030405403171.csv", 
                   skip = 13) %>% 
  rename(co2_ppm = `CO2 (Dis)@110118B`, q_ls = `Discharge@110118B`,
         temp_water_c = `Water Temp@110118B`, datetime = TimeStamp) %>% 
  mutate(year= year(datetime), month= month(datetime), day= day(datetime), hour= hour(datetime),
         date_hour = make_datetime(year, month, day, hour, 0, 0)) %>% 
  select(date_hour, co2_ppm, q_ls, temp_water_c) %>% 
  summarise(across(everything(), mean), .by ="date_hour") %>% 
  mutate(site= "k18")


c4_co2 <- read_csv("empirical data/river data/streampulse/_raw/krycklan/raw/C4_Q_TW_CO2_2017-2022.Time_Series_Data.2024030405272149.csv", 
                    skip = 13) %>% 
  rename(co2_ppm = `CO2 (Dis)@110104_C4`, q_ls = `Discharge@110104_C4`,
         temp_water_c = `Water Temp@110104_C4`, datetime = TimeStamp) %>% 
  mutate(year= year(datetime), month= month(datetime), day= day(datetime), hour= hour(datetime),
         date_hour = make_datetime(year, month, day, hour, 0, 0)) %>% 
  select(date_hour, co2_ppm, q_ls, temp_water_c) %>% 
  summarise(across(everything(), mean), .by ="date_hour") %>% 
  mutate(site= "k4")

c6_co2 <- read_csv("empirical data/river data/streampulse/_raw/krycklan/raw/C6_Q_TW_CO2.Time_Series_Data.2024030405322643.csv", 
                   skip = 13) %>% 
  filter(as.Date(TimeStamp) < as.Date("2021-05-01")) %>% 
  rename(co2_ppm = `CO2 (Dis)@110106_C6`, q_ls = `Discharge@110106_C6`,
         temp_water_c = `Water Temp@110106_C6`, datetime = TimeStamp ) %>% 
  mutate(year= year(datetime), month= month(datetime), day= day(datetime), hour= hour(datetime),
         date_hour = make_datetime(year, month, day, hour, 0, 0)) %>% 
  select(date_hour, co2_ppm, q_ls, temp_water_c) %>% 
summarise(across(everything(), mean), .by ="date_hour") %>% 
  mutate(site= "k6")

#turn the do dataset into hourly
krycklan_do_hourly <- krycklan_do %>% 
  mutate(year= year(datetime_local_ok), month= month(datetime_local_ok), day= day(datetime_local_ok), hour= hour(datetime_local_ok),
         date_hour = make_datetime(year, month, day, hour, 0, 0)) %>% 
  summarise(across(c(do_con, temp, do_sat), mean), .by= c("date_hour", "site")) %>% 
  mutate(site= str_remove(site, "_up")) 

#join with co2 and o2"
co2_o2_krycklan <- bind_rows(c4_co2, c6_co2, c7_co2) %>% 
  left_join(krycklan_do_hourly, by= c("date_hour", "site")) %>% 
  mutate(siteID = site,
         dateTimeUTC = date_hour,
         SpecCond_uScm = NA,
         WaterTemp_C = temp_water_c,
         Depth_m = NA,
         Light_lux = NA,
         Discharge_m3s = q_ls/1000,
         AirPressCombined_atm = NA,
         CO2_mmolm3 = ppm_to_mmolm3(co2_ppm, 1, WaterTemp_C),
         O2_mmolm3 = do_con/32*1000,
         CO2eq_mmolm3 = ppm_to_mmolm3(air_CO2, 1, WaterTemp_C),
         O2eq_mmolm3 = osat_grace(WaterTemp_C, 1, sal =0.001)/32*1000,
         CO2dep_mmolm3 = CO2_mmolm3 - CO2eq_mmolm3,
         O2dep_mmolm3 = O2_mmolm3 - O2eq_mmolm3,
         CDOM_ppb=NA,
         pH=NA, 
         Turbidity_NTU=NA, 
         Turbidity_FNU=NA, 
         SpecCond_uScm=NA, 
         Nitrate_mgL=NA ) %>% 
  select(siteID, dateTimeUTC, WaterTemp_C, Depth_m,  Light_lux, Discharge_m3s, AirPressCombined_atm, CO2_mmolm3:O2dep_mmolm3, 
         CDOM_ppb, pH, Turbidity_NTU, Turbidity_FNU, SpecCond_uScm, Nitrate_mgL) %>% 
  filter(as.Date(dateTimeUTC) < as.Date("2018-07-01"))



#quick look  
co2_o2_krycklan %>% 
  drop_na(CO2dep_mmolm3, O2dep_mmolm3) %>% 
  mutate(year= year(dateTimeUTC)  ) %>%
  arrange(desc(year)) %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3, color=month(dateTimeUTC) ))+
  geom_point(alpha=.1)+
  geom_abline(slope=-1, intercept = 0, linetype= 2)+ 
  facet_wrap(~siteID, scales = "free")+
  theme_classic()+
  theme(legend.position = "top")+
  labs(color= "Year")+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


  
#put data in the right format
sp_data_out <- sp_data  %>% 
  left_join(sp_site_data %>% select(siteID, AirPresElev_atm, latitude, longitude), by = "siteID") %>% 
  mutate(SpecCond_uScm = if_else(is.na(SpecCond_uScm) == TRUE, SpecCond_mScm/1000, SpecCond_uScm),
         AirPressCombined_atm = if_else(is.na(AirPres_atm) == TRUE, AirPresElev_atm, AirPres_atm),
         CO2_mmolm3 = ppm_to_mmolm3(CO2_ppm, AirPressCombined_atm, WaterTemp_C),
         O2_mmolm3 = DO_mgL/32*1000,
         CO2eq_mmolm3 = ppm_to_mmolm3(air_CO2, AirPressCombined_atm, WaterTemp_C),
         O2eq_mmolm3 = osat_grace(WaterTemp_C, AirPressCombined_atm, sal =0.001)/32*1000,
         CO2dep_mmolm3 = CO2_mmolm3 - CO2eq_mmolm3,
         O2dep_mmolm3 = O2_mmolm3 - O2eq_mmolm3 ) 


sp_data_out %>% 
  select(c(AirPressCombined_atm, WaterTemp_C, CO2_ppm )) %>%  
  summarise_all(funs(sum(is.na(.))))

sp_data_out %>% 
  select(c(AirPressCombined_atm, WaterTemp_C, CO2_ppm )) %>%  #
  summarise_all(funs(max(.)))




sp_data_cleaned <- sp_data_out %>% 
  mutate(remove = case_when(siteID == "QS" & WaterTemp_C < 2 ~ "yes",
                            siteID == "NHC" & CO2_ppm > 15000 ~ "yes",
                            siteID == "NHC" & DO_mgL > 20 ~ "yes",
                            siteID == "Drain" & CO2_ppm < 2000~ "yes",
                            siteID == "Drain" & O2dep_mmolm3 > -20 ~ "yes",
                            siteID == "BEC" & dateTimeUTC > as.Date("2017-10-01") ~ "yes",
                            siteID == "BEC" & CO2_ppm > 10000 ~ "yes",
                            siteID == "Eno" ~ "yes",
                            siteID == "BEC" & dateTimeUTC > as.Date("2018-02-01") ~ "yes",
                            siteID == "Drain" & dateTimeUTC < as.Date("2019-01-01") ~ "yes",
                            siteID == "Drain" & dateTimeUTC < as.Date("2019-07-01") & CO2_ppm < 3000 ~ "yes",
                            siteID == "NR1000" & dateTimeUTC > as.Date("2020-07-01") ~ "yes",
                            siteID == "NR1000" & CO2dep_mmolm3 > 280 ~ "yes",
                            siteID == "SF2500" & dateTimeUTC > as.Date("2018-06-01") & dateTimeUTC < as.Date("2018-07-23") ~ "yes",
                            siteID == "SF2500" & dateTimeUTC > as.Date("2019-08-03") & dateTimeUTC < as.Date("2019-08-13") ~"yes",
                            siteID == "SF2800" & dateTimeUTC > as.Date("2019-09-13") ~ "yes",
                            siteID == "SF2800" & dateTimeUTC > as.Date("2018-10-21") & dateTimeUTC < as.Date("2018-11-23") ~ "yes",
                            TRUE ~ "no") ) %>% 
  filter(remove == "no") %>% 
  select(siteID, dateTimeUTC, WaterTemp_C, Depth_m,  Light_lux, Discharge_m3s, AirPressCombined_atm, CO2_mmolm3:O2dep_mmolm3, 
         CDOM_ppb, pH, Turbidity_NTU, Turbidity_FNU, SpecCond_uScm, Nitrate_mgL)

##### check some data 
sp_data_cleaned %>% 
  filter(siteID == "BRW") %>% 
  ggplot(aes(CO2dep_mmolm3, O2dep_mmolm3, color=WaterTemp_C))+
  geom_point(alpha=.2)+
  scale_color_viridis_c()

##### check some data 
co2_o2_krycklan %>% 
  filter(siteID == "k7") %>% 
  drop_na(O2dep_mmolm3) %>% 
  ggplot(aes(dateTimeUTC, CO2dep_mmolm3, color=WaterTemp_C))+
  geom_point(alpha=.2)+
  scale_color_viridis_c()


#export the time series dataset
bind_rows(sp_data_cleaned, miellajokka, co2_o2_krycklan)  %>% 
  write_csv("empirical data/river data/streampulse/stream_dataset.csv")


## get USGS data for each site
library(dataRetrieval)


sites_usgs <- sp_site_data %>% drop_na(USGSgageID)

for(i in seq_along(sites_usgs$USGSgageID)){
  
  dat_sites <- readWQPdata(
    siteid = paste0("USGS-",sites_usgs$USGSgageID[i]),
    dataProfile = "resultPhysChem")
  
  
  write_csv(dat_sites, paste0("data/USGS_data/chemistry_usgs-",sites_usgs$USGSgageID[i], ".csv"))
  
  print(paste("done", i))
  
}



for(i in seq_along(sites_usgs$USGSgageID)){
  
  dat_sites <- readNWISuv(siteNumbers = as.character(sites_usgs$USGSgageID[i]),
                          parameterCd = "00060",
                          startDate = "2010-01-01",
                          endDate = "2021-12-31")
  
  
  write_csv(dat_sites, paste0("data/USGS_data/discharge_usgs-",sp_site_data$USGSgageID[i], ".csv"))
  
  print(paste("done", i))
  
}


chemistry_all <- list.files(path = "data/USGS_data/chemistry/",
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

chemistry_all %>% 
  write_csv("data/stream_datasets/chem_usgs.csv")


discharge_all <- list.files(path = "data/USGS_data/discharge/",
                            pattern="*.csv", 
                            full.names = T) %>% 
  map_df(~read_csv(., col_types = cols(.default = col_character()))) %>% 
  mutate(date_time= as.character(dateTime),
         site= as.character(site_no),
         discharge = as.character(X_00060_00000),
         .keep = "none")

discharge_all %>% 
  write_csv("data/stream_datasets/discharge_usgs.csv")
