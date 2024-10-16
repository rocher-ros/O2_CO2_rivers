
#make an artificial river, with radiation based on longitude and latitude, and a temperature curve that goes up at day


# OSAT FUNCTION 
#oxygen saturation function based on temp, salinity and pressure as in Grace et al., 2015 L&O
#pressure in atmospheres and temp in celsius
o2_sat <- function(temp, bp, sal){ 
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
  
  return(sato)
}

#co2 in equilibrium with the atmosphere, depending on temperature (celsius) and pressure (atm), as mmol m-3
#original source of constants:
#Millero, F. (2010). Carbonate constants for estuarine waters Marine and Freshwater
#Research 61(2), 139. As amended by Orr et al. 2015.
co2_sat <- function(temp, bp, sal){
  
  pCO2.uatm = 400 #uatm, change in a few years :S
  
  Bar.pressure = bp*101325 #change to Pa
  
  Kh2 = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp)) + 23.3585 * log( (273.15 + temp) / 100 ) + 
              sal*(0.023517 - 0.023656*(273.15 + temp)/100 + 0.0047036*((273.15 + temp)/100)^2)) #mol L-1 Pa-1
  
  pCO2.uatm*Kh2/Bar.pressure*101.325*1e+3 #mmol m-3
  
}




### K FUNCTION ###
# calculate KO2 from K600, corrected for temperature (units are day-1)
# KO2 = schmidt.O2/schmidt.600 * K600 (Wanninkhof 1992 JGR)
K_O2 <- function (temp,K600) {
  schmidt.O2 <-1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)
  
  K600/(600/(schmidt.O2))^(0.5)
}

#calculate K_CO2 using tame and a given K600
K_CO2 <- function (temp,K600) {
  schmidt.CO2 <- (1911.1-(118.11*(temp))+(3.4527*(temp^2))-(0.04132*(temp^3)))#schmidt number of the water, depending on temperature
  
  K600/(600/(schmidt.CO2))^(0.5)
}

## # # # #  function to calculate co2 from alk and DIC USING PHREEQC 
# using an oldscript form Gerard rocher-Ros
# some changes taken from: https://doi.org/10.1016/j.mex.2021.101430
dic_equil_co2 <- function(dic, alk, tc){
  
  #Give input data to PHREEQC
  uni <- "umol/l" #units of C
  
  #Make strings with data
  unis <- paste("units ",uni,sep=" ")
  temps <- paste("temp ",as.character(tc),sep=" ")
  alks <- paste("Alkalinity ",as.character(alk),sep=" ")
  dics = paste("C(4) ",as.character(dic),sep=" ")
  
  
  #make input string for phreeqc
  is <- c("SOLUTION 1 water",
          unis,
          temps,
          dics,
          alks,
          "SELECTED_OUTPUT",
          "-file spec.sel",
          "-totals C(4) C(-4)",
          "CO2(g)",
          "-molalities CO2 CO3-2 HCO3-"
  )
  
  #convert indata to data frame
  isd <- as.data.frame(is)
  
  #load the phreeqc.dat database
  phrLoadDatabaseString(phreeqc.dat)
  
  #run phreeqc
  phrRunString(is)
  
  #retrieve selected_output as a list of data.frame
  os <- phrGetSelectedOutput()
  
  #convert results to data frame
  osd <- as.data.frame(os$n1)
  
  data_out <- tibble(  #get  the calculated  data from Phreeqc, bring back to used units
    co2.mmol.m3 = osd$m_CO2.mol.kgw.*1e+6,
    hco3.mmol.m3 = osd$m_HCO3..mol.kgw.*1e+6,
    co3.mmol.m3 = osd$m_CO3.2.mol.kgw.*1e+6,
    pH = osd$pH) 
  
  return(tibble(co2.mmolm3 = data_out$co2.mmol.m3,
                hco3.mmolm3= data_out$hco3.mmol.m3,
                co3.mmolm3= data_out$co3.mmol.m3,
                pH= data_out$pH))
}


#radiation function 
hourly_radiation <- function (jday, solrad, hour, lon, lat){ #from package TrenchR
  rd <- 180/pi
  RevAng <- 0.21631 + 2 * atan(0.967 * tan(0.0086 * (-186 + jday)))
  DecAng <- asin(0.39795 * cos(RevAng))
  f <- (279.575 + 0.9856 * jday)/rd
  ET <- (-104.7 * sin(f) + 596.2 * sin(2 * f) + 4.3 * sin(3 * f) - 12.7 * sin(4 * f) - 429.3 * cos(f) - 2 * cos(2 *f) + 19.3 * cos(3 * f))/3600
  LC <- 1/15 * (15 - lon %% 15)
  hour_sol <- 12 + LC + ET
  W <- pi * (hour - hour_sol)/12
  Ws <- acos(-tan(lat/rd) * tan(DecAng))
  d <- 0.409 + 0.5016 * sin(Ws - 1.047)
  b <- 0.6609 - 0.4767 * sin(Ws - 1.047)
  rG <- pi/24 * (d + b * cos(W)) * (cos(W) - cos(Ws))/(sin(Ws) -  Ws * cos(Ws))
  
  rG[rG < 0] = 0
  
  return(rG * solrad)
  
}

# temperature function
hourly_temp <- function(T_max, T_min, hour){ #from TrenchR package
  #stopifnot(t >= 0, t <= 24, T_max >= T_min)
  W <- pi/12
  gamma <- 0.44 - 0.46 * sin(0.9 + W * hour) + 0.11 * sin(0.9 +  2 * W * hour)
  T_max * gamma + T_min * (1 - gamma)
  
}

#estimate atompsheric pressure (in atm) from elevation (m)
estimate_pressure <- function(elevation){
  pressure <- (1-(.0000225577*elevation))^5.25588
  return(pressure)
}


create_river <- function(jday, lat, lon, elev, radiation.day, temp.mean, temp.shift, res= dt){
  
#estimate atm.pressure based on elevation, 
atm.press= estimate_pressure(elev) #atm
temp.var = temp.shift/2


#First we start with an empty time series vector
ts <- seq(as.POSIXct("2020-07-01", tz = "UTC"), 
          as.POSIXct("2020-07-15", tz = "UTC"), 
          by = paste(res, "min"))

#remove the last value because it starts a new day
ts <- ts[-length(ts)]

df_init <- tibble(date.time= ts,
                         hour= hour(date.time) + minute(date.time)/60,
                         date= as.Date(date.time)) %>% 
  mutate(radiation = hourly_radiation(jday, radiation.day, hour, lon, lat),
         temp.air = hourly_temp(T_max = temp.mean + temp.var, T_min= temp.mean - temp.var, hour= hour),
         temp.water = hourly_temp(T_max = temp.mean + temp.var/2, T_min= temp.mean - temp.var/2, hour= hour-2),
         atm.press = atm.press)  #current approach with water temp. a bit buffered from air and a 2h lag, just eye balling it for now

return(df_init)
}


plots_function <- function(data, column){
  
  plot_o2 <- ggplot(data)+
    geom_point(aes(date.time, o2.mod/o2.air*100, color= .data[[column]]))+
    scale_color_viridis_c()+
    labs(y="O2 (%sat)")+
    theme_bw()
  
  plot_dic <- ggplot(data)+
    geom_point(aes(date.time, dic.mod, color= .data[[column]]))+
    scale_color_viridis_c()+
    labs(y="DIC (mmol m3)")+
    theme_bw()
  
  plot_ph <- ggplot(data)+
    geom_point(aes(date.time, pH, color= .data[[column]]))+
    scale_color_viridis_c()+
    theme_bw()
  
  plot_co2 <- ggplot(data)+
    geom_point(aes(date.time, co2.mod/co2.air*100, color= .data[[column]]))+
    scale_color_viridis_c()+
    labs(y="CO2 (%sat)")+
    theme_bw()
  
  
  #CO2:O2 plot
  plot_elipse_CO2 <- data %>% 
    #arrange(desc(.data[[column]])) %>% 
    ggplot(aes(co2.mod-co2.air, o2.mod-o2.air, color= .data[[column]]))+
    geom_path(linewidth=2, alpha=.3)+
    geom_abline(slope=-1, intercept = 0, linetype= 2)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    scale_color_viridis_c()+
    theme_bw()+
    labs(x="departure CO2", y= "departure O2")
  
  #DIC:O2 plot
  plot_elipse_DIC <- data %>% 
    #arrange(desc(.data[[column]])) %>% 
    ggplot(aes(dic.mod-co2.air, o2.mod-o2.air, color= .data[[column]]))+
    geom_path(linewidth=2, alpha=.3)+
    geom_abline(slope=-1, intercept = 0, linetype= 2)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    scale_color_viridis_c()+
    theme_bw()+
    labs(x="departure DIC", y= "departure O2")
  
  plot_elipse <- plot_elipse_CO2 + plot_elipse_DIC +
    plot_layout(guides = "collect") & theme(legend.position = "top")
  
  ggsave(filename = paste0("plots/model_output/", column, "_elipse.png"), plot = plot_elipse, scale = .8)
  
  
  variables_plot <- plot_o2+ plot_dic + plot_ph + plot_co2+ plot_layout(guides = "collect") & theme(legend.position = "top")
  
  ggsave(filename = paste0("plots/model_output/", column, "_variables.png"), plot = variables_plot)
  
  return(plot_elipse)
  
}

