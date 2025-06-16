# Model to reproduce dissolved inorganic carbon system and O2 under the effect of multiple drivers in streams.
# Model developed by G.Rocher-Ros, based on previous models by Stets et al (2017) and Tank & hall (2014)
# This script includes the main model function. There is another script with auxiliary functions needed to call it

source("scripts/model_o2_co2/auxiliary_functions.R")


# Install and Load libraries ----

# List of all packages needed
package_list <- c('tidyverse', 'phreeqc')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

# Now load all the packages
lapply(package_list, require, character.only = TRUE)

#for debugging
params <- tibble(GPP.day= 3/32*1e+3, ER.day = 8/32*1e+3, PQ= 1.2, RQ= 1.2,  bigK = 20, alk.mmolm3= 500, q.m3s= 1.3, 
                depth.m = .5, gw.frac= 0.05, gw.O2conc = 150, gw.CO2conc = 550 ) 


# Model call main function -----
co2_o2_sim <- function( params ){
  
  #temporal resolution of the model
  dt= 10 # minutes
  
  #pull out all parameters into values for better control
  GPP.day= params$GPP.day
  ER.day = params$ER.day
  PQ = params$PQ
  RQ = params$RQ
  bigK= params$bigK
  q.m3s= params$q.m3s
  depth.m = params$depth.m
  slope = params$slope
  alk.mmolm3= params$alk.mmolm3
  gw.frac= params$gw.frac
  gw.O2conc = params$gw.O2conc
  gw.DICconc= params$gw.CO2conc+alk.mmolm3
  
  
  #create a time series with temp and radiation day-night changes, based on location and some set parameters
  df_init <- create_river(
    res= dt, #resolution of time series
    jday= 110, #julian day
    lat=42,# degrees
    lon= 2, #degrees
    elev= 100, #meters
    radiation.day = 2000, # W m2
    temp.mean = 20, # degrees celsius, mean temp
    temp.shift = 10) #day/night shift in air temp
  
  #in case you don't provide depth, can be estimated from Q
  if(is.null(depth.m) == TRUE) {
    depth.m = exp(-0.895 + 0.294 * log(q.m3s)) #from raymond et al 2012
   }  
  
  
  #in case you don't have K, can be estimated from hydraulics
  if(is.null(bigK) == TRUE) {
   
    vel.ms = exp(0.12 * log(q.m3s) - 1.06) #using USGS Q-V relationship,
    width.m = exp(2.1 + 0.47 * log(q.m3s)) #from Liu et al 2022 
    k600.md = 2841 * slope * vel.ms + 2.02 # m d-1 #Raymond et al 2012
    K600.d = k600.md/ depth.m # d-1
    message("if you don't have a K you need to provide Q and slope")
  }  


  #initial daily values. first calculate total radiation and then atatch the parameters and other stuff
  df_init_daily <- df_init %>% 
    summarise(radiation.day = sum(radiation), .by= "date") %>% 
    mutate(alk.mmolm3 = alk.mmolm3,
           gw.O2conc = gw.O2conc, 
           gw.DICconc = gw.DICconc,
           PQ = PQ,
           RQ = RQ,
           q.m3s = q.m3s,
           gw.frac= gw.frac,
           depth.m= depth.m,#exp(-0.895 + 0.294 * log(q.m3s)), #from raymond et al 2012
           vel.ms = exp(0.12 * log(q.m3s) - 1.06), #using USGS Q-V relationship,
           width.m = exp(2.1 + 0.47 * log(q.m3s)), #from Liu et al 2022 
           K600.d = bigK, # d-1
           k600.md = K600.d*depth.m, 
           footprint.m = -log(1 - 0.8)*vel.ms*3600*24/K600.d, #80% gas footprint
           area.m2 = footprint.m*width.m, #m2
           GPP.day= GPP.day,
           ER.day = ER.day)
  

 #dataframe at the dt resolution, currently set at 10min, with  flux values in min
  df <- df_init %>% 
    left_join(df_init_daily, by= "date") %>% 
    mutate(q.m3min = q.m3s*60, #m3 min-1
           q.gw = q.m3min*gw.frac,
           K.O2 = K_O2(temp.water, K600.d)/24/60, # from d-1 to min-1
           o2.air = o2_sat(temp.water, atm.press, sal= 0)/32*1e+3, #function is in mgO2/L, convert to mmol/m3
           K.CO2 = K_CO2(temp.water, K600.d)/24/60, #from d-1 min-1
           co2.air = co2_sat(temp.water, atm.press, sal= 0),
           GPP = GPP.day*(radiation/radiation.day)/dt, # mmol O2 m-2 min-1
           ER = ER.day/24/60, # mmol O2 m-2 min-1
           o2.mod = NA,   #create empty columns for values to fill
           dic.mod = NA,
           co2.mod = NA,
           pH = NA) 
  

  #set some starting conditions, atmospheric equilibrium for everything
  df$o2.mod[1] = df$o2.air[1]
  df$dic.mod[1] = alk.mmolm3 + df$co2.air[1]
  df$co2.mod[1] = df$co2.air[1]
  
  
  #Main loop to run the simulation
  for(i in 2:nrow(df)){
    
    #Main Eq for oxygen, following Hall & Tank (2005)
    df$o2.mod[i] = df$o2.mod[i-1] +  #O2 from previous time step
      (df$q.gw[i]*gw.O2conc)/(df$area.m2[i]*df$depth.m[i])*dt - (df$q.gw[i]*df$o2.mod[i-1])/(df$area.m2[i]*df$depth.m[i])*dt + #groundwater inputs
      df$GPP[i]*dt/df$depth.m[i] +#input from GPP
      df$ER[i]*dt/df$depth.m[i] -  #output to ER
      df$K.O2[i]*(df$o2.mod[i-1] - df$o2.air[i])*dt #exchange with the atmosphere
    
    #modelled DIC, changing from metabolism and GW inputs (now disabled)
    df$dic.mod[i] = df$dic.mod[i-1] +  #dic from previous time step
      (df$q.gw[i]*gw.DICconc)/(df$area.m2[i]*df$depth.m[i])*dt - (df$q.gw[i]*df$dic.mod[i-1])/(df$area.m2[i]*df$depth.m[i])*dt - #groundwater inputs from Hall & Tank 2005
      df$GPP[i]*1/df$PQ[i]*dt/df$depth.m[i] -# output to GPP
      df$ER[i]*df$RQ[i]*dt/df$depth.m[i] -  # input from ER, is negative as ER is negative, so is plus
      df$K.CO2[i]*(df$co2.mod[i-1] - df$co2.air[i])*dt #exchange with the atmosphere, with co2 from previous ts

    
    #calculate the DIC fractions with the given DIC and fixed alk
    dic.equil <- dic_equil_co2(dic= df$dic.mod[i], alk = alk.mmolm3, tc = df$temp.water[i])
    

    
    #Update the modelled co2 from the DIC equlibration
    df$co2.mod[i] = dic.equil$co2.mmolm3   #CO2 from equilibrated  dic
    

    
    #save the calculated pH
    df$pH[i] = dic.equil$pH
    
  }
  
  df_out <- df %>% 
    slice(-c(1:300))
  
  return(df_out)
}





