
# LOAD PACKAGES AND READ FILES ----
## Install and Load libraries ----

# List of all packages needed
package_list <- c('tidyverse', 'patchwork', 'sf', 'lmodel2', 'lwgeom')

# Check if there are any packacges missing
packages_missing <- setdiff(package_list, rownames(installed.packages()))

# If we find a package missing, install them
if(length(packages_missing) >= 1) install.packages(packages_missing) 

# Now load all the packages
lapply(package_list, require, character.only = TRUE)



#Load all custom functions for the model
source("scripts/model_o2_co2/model_dic_o2.R")


#functions to calculate ellipse and regression metrics ----
# Modified from Vachon et al., 2020 GCB
# To  calculate the stats we need a minimum amount of days, so to avoid incomplete periods
# here I have set a minimum of 48 observations which correspond to 2 days, 
elipse <- function(CO2, O2){
  
  if(length(CO2) > 48){
    corMat = cor(cbind(CO2, O2))
    covMat <- var(cbind(CO2, O2))
    evals <- eigen(covMat)$values
    ell.len <- 2*sqrt(5.991*evals)
    
    out <- tibble(width = ell.len[2], #output is the length and width of the ellipse
                  length = ell.len[1]) 
  }else{
    out <- tibble(width = NA, length = NA) 
  }
  out
}

#Type II regression for slope calculation 
regression_metrics <- function(CO2, O2){
  
  if(length(CO2) > 48){  
    reg <- lmodel2(O2 ~ CO2, nperm = 99)
    slope <- reg$regression.results[2,3]
    corr <- reg$r
    
    out <- tibble(EQ= 1/abs(slope), #the output is the EQ (inverse of slope)
                  corr= corr) #and the pearson r
  }else{
    out <- tibble(EQ = NA, corr = NA) 
  } 
  out
  
}


#load metab data from Appling
daily_metab <- read_delim("prepared data/river data/Appling2019/daily_predictions.tsv") %>% 
  filter(GPP>=0, ER <= 0, GPP.Rhat < 1.1, ER.Rhat < 1.1, K600.Rhat < 1.1)

site_info <- read_delim("prepared data/river data/Appling2019/site_data.tsv")

site_catchments <- read_sf("prepared data/river data/Appling2019/catchment_shapefile/catchment_shapefile.shp")
outlet_catchments <- read_sf("prepared data/river data/Appling2019/points_shapefile/points_shapefile.shp")

#read files with catchment network properties downloaded before
catchment_lengths <- read_csv("prepared data/river data/USGS_data/catchment_lengths_properties.csv") %>% 
  mutate(avg_discharge_m3s= 0.028316832*avg_dicharge) %>%  #convert to m3s from cfs
  select(-avg_dicharge)

names(catchment_lengths)

sf_use_s2(FALSE)

# INTIAL DATA WRANGLING ----

#make a new tibble with catchment areas from the shp file from Appling 2019
catch_areas <- tibble(site_name = site_catchments$site_nm) %>% 
  mutate( catchment_area = st_area(site_catchments) %>% as.numeric())

#perform summary stats for each site with the daily metabolism file
#also join other useful files, catchment areas and network lengths, and  estimate the %gw inputs in each site
site_metab <- daily_metab %>% 
  summarise(across(c(GPP, ER, K600, DO.obs, DO.amp, depth, discharge, temp.water, shortwave, DO.tdist80, velocity), 
                   list(mean= ~mean(., na.rm= TRUE), 
                        median= ~median(., na.rm= TRUE), 
                        cv= ~sd(., na.rm= TRUE)/mean(., na.rm= TRUE)*100)), .by = "site_name") %>% 
  left_join(catch_areas, by = "site_name") %>% 
  left_join(catchment_lengths, by= c("site_name"="site_nm")) %>% 
  mutate(groundwater_m3s = avg_discharge_m3s- (tot_length-DO.tdist80_median/1000)/tot_length*avg_discharge_m3s,
         gw_frac= groundwater_m3s/avg_discharge_m3s*100,
         spQ_mmday = discharge_mean*3600*24/catchment_area*1000) 

write_csv(site_metab, "prepared data/river data/Appling2019/site_avgs_gwinputs.tsv")

# INITIAL VISUALIZATION OF THE DATA ----

site_metab %>% 
  filter(tot_area/(catchment_area/1e+6) > .5, tot_area/(catchment_area/1e+6) < 5,
         gw_frac <= 100, gw_frac > .01,
         discharge_median> 0.07, discharge_mean < 1e+4) %>% 
  ggplot(aes(discharge_median, gw_frac))+
  geom_smooth(color="gray9")+
  geom_point(aes(color=spQ_mmday))+
  scale_x_log10(labels = scales::number, breaks= c(0.1, 1, 10, 100, 1000, 10000))+
  scale_y_log10(labels = scales::number)+
  scale_color_viridis_c(trans= "log10")+
  theme_bw()+
  labs(x= expression(Discharge~(m^3~s)), y = "% discharge from groundwater within the 80% footprint reach",
       color= "Specific discharge (mm d-1)")+
  theme(legend.position = "top")




# RUN THE MODEL ACROSS THE RCC ----
#prepare data for the model, here we exclude extreme values for GPP, ER, Q, as well as sites where the
# GIS extracted catchment area doesn't match the one for the datasets
data_for_model <- site_metab %>% 
  filter(GPP_mean < 10, ER_mean > -20, discharge_mean < 500, discharge_mean > 0.006,
         tot_area/(catchment_area/1e+6) > .5, tot_area/(catchment_area/1e+6) < 5,
         gw_frac <= 100) %>% 
  mutate(NEP_mean= ER_mean + GPP_mean,
         NEP_median = ER_median + GPP_median) %>% 
  left_join(site_info, by = "site_name")


#extract trends of parameters with discharge
lw_GPP <- loess(GPP_mean ~ log(discharge_mean), data = data_for_model)
lw_ER <- loess(ER_mean ~ log(discharge_mean), data = data_for_model)
lw_NEP <- loess(NEP_mean ~ log(discharge_mean), data = data_for_model)
lw_K <- loess(K600_mean ~ log(discharge_mean), data = data_for_model)
lw_depth <- loess(depth_mean ~ log(discharge_mean), data = data_for_model)
lw_gw <- loess(gw_frac ~ log(discharge_mean), data = data_for_model)

# predict fitted values for each value and each parameter
fit_GPP <- data.frame(predict(lw_GPP, se = TRUE))
fit_ER <- data.frame(predict(lw_ER, se = TRUE))
fit_NEP <- data.frame(predict(lw_NEP, se = TRUE))
fit_K <- data.frame(predict(lw_K, se = TRUE))
fit_depth <- data.frame(predict(lw_depth, se = TRUE))
fit_gw <- data.frame(predict(lw_gw, se = TRUE))


# make a new tibble with all variables across dicharge
RCC_smoothed <- tibble(
             discharge = data_for_model$discharge_mean,
             depth.m = data_for_model$depth_mean,
             GPP = data_for_model$GPP_mean,
             GPP.fit = fit_GPP$fit,
             GPP.upperBound = fit_GPP$fit + 2 * fit_GPP$se.fit,
             GPP.lowerBound = fit_GPP$fit - 2 * fit_GPP$se.fit,
             ER = data_for_model$ER_mean,
             ER.fit = fit_ER$fit,
             ER.upperBound = fit_ER$fit + 2 * fit_ER$se.fit,
             ER.lowerBound = fit_ER$fit - 2 * fit_ER$se.fit,
             NEP = data_for_model$NEP_mean,
             NEP.fit = fit_NEP$fit,
             NEP.upperBound = fit_NEP$fit + 2 * fit_NEP$se.fit,
             NEP.lowerBound = fit_NEP$fit - 2 * fit_NEP$se.fit,
             K600 = data_for_model$K600_mean,
             K600.fit = fit_K$fit,
             K600.upperBound = fit_K$fit + 2 * fit_K$se.fit,
             K600.lowerBound = fit_K$fit - 2 * fit_K$se.fit,
             gw = data_for_model$gw_frac,
             gw.fit = fit_gw$fit,
             gw.upperBound = fit_gw$fit + 2 * fit_gw$se.fit,
             gw.lowerBound = fit_gw$fit - 2 * fit_gw$se.fit)

plot_nep <- ggplot(RCC_smoothed)+
  geom_ribbon(aes( discharge, ymin = NEP.lowerBound, ymax= NEP.upperBound), alpha=.3, fill= "yellow4")+
  geom_point(aes( discharge, NEP), alpha=.5, color="yellow3")+
  geom_line(aes( discharge, NEP.fit), color="yellow4")+
  geom_hline(yintercept = 0)+
  scale_x_log10(labels = scales::number)+
  theme_classic()+
  labs(x= "", y= expression(NEP~(g~O[2]~m^-2~d^-1)))

plot_gpp <- ggplot(RCC_smoothed)+
  geom_ribbon(aes( discharge, ymin = GPP.lowerBound, ymax= GPP.upperBound), alpha=.3, fill= "chartreuse3")+
  geom_point(aes( discharge, GPP), alpha=.5, color="chartreuse4")+
  geom_line(aes( discharge, GPP.fit), color="chartreuse4")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_log10(labels = scales::number)+
  theme_classic()+
  labs(x= "", y=expression(GPP~(g~O[2]~m^-2~d^-1)))


plot_er <- ggplot(RCC_smoothed)+
  geom_ribbon(aes( discharge, ymin = ER.lowerBound, ymax= ER.upperBound), alpha=.3, fill= "red4")+
  geom_point(aes( discharge, ER), alpha=.5, color="red4")+
  geom_line(aes( discharge, ER.fit), color="red3")+
  geom_hline(yintercept = 0)+
  scale_x_log10(labels = scales::number)+
  theme_classic()+
  labs(x= "", y= expression(ER~(g~O[2]~m^-2~d^-1)))


plot_k <- RCC_smoothed %>% 
  filter(K600<60) %>% 
ggplot()+
  geom_ribbon(aes( discharge, ymin = K600.lowerBound, ymax= K600.upperBound), alpha=.3, fill= "gray60")+
  geom_point(aes( discharge, K600), alpha=.5, color="gray30")+
  geom_line(aes( discharge, K600.fit), color="gray30")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_log10(labels = scales::number)+
  theme_classic()+
  labs(x= "", y= expression(K[600]~(d^-1)))

plot_gw <-
  RCC_smoothed %>% 
  ggplot()+
  geom_ribbon(aes( discharge, ymin = gw.lowerBound, ymax= gw.upperBound), alpha=.3, fill= "dodgerblue2")+
  geom_point(aes( discharge, gw), alpha=.5, color="dodgerblue4")+
  geom_line(aes( discharge, gw.fit), color="dodgerblue4")+
  scale_x_log10(labels = scales::number)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x= expression(Discharge~(m^3~s^-1)), y= "% Groundwater")

plot_vars_rcc <- plot_gpp + plot_er + plot_k + plot_gw + 
  plot_layout(ncol =1 ) +
  plot_annotation(tag_levels = 'a')

plot_vars_rcc

ggsave("plots/main/fig4_vars_rcc.png", plot = plot_vars_rcc,  width= 6, height = 8)

params_rcc <- tibble(
  discharge_mean= 10^(seq(from = log10(min(data_for_model$discharge_mean)), 
                          to = log10(max(data_for_model$discharge_mean)), length.out=50)) ) %>% 
  mutate(GPP.day = predict(lw_GPP, newdata =. )/32*1e+3,
         ER.day = predict(lw_ER, newdata=.)/32*1e+3,
         bigK = predict(lw_K, newdata=.),
         depth.m = predict(lw_depth, newdata=.),
         gw.frac = predict(lw_gw, newdata=.)/100,
         PQ = 1.2,
         RQ = 1.2,
         alk.mmolm3 = 2221,
         gw.O2conc = 5.54/32*1e+3, 
         gw.CO2conc = 772) %>% 
  rename(q.m3s= discharge_mean)





## run the simulation with the set parameters
sim_out <- params_rcc %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim, .id= "id")


plot_ellipse_rcc <- sim_out %>% 
  mutate(q.categories = case_when(q.m3s < 0.1 ~ "very small",
                                  q.m3s >= 0.1 & q.m3s < 1 ~ "small",
                                  q.m3s >= 1 & q.m3s < 10 ~ "medium",
                                  q.m3s >= 10 & q.m3s < 100 ~ "big",
                                  q.m3s >= 100 ~ "very big") %>% 
           as_factor() %>% 
           fct_reorder(q.m3s, mean, .desc = TRUE)) %>% 
  filter(date > as.Date("2020-07-13")) %>% 
  ggplot(aes(co2.mod- co2.air, o2.mod-o2.air, color= q.m3s))+
  geom_point(size= .4, alpha=.4)+
  geom_abline(slope=-1, intercept = 0, linetype= 2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_y_continuous(breaks= c(-50, -25, 0, 25), limits = c(-60, 25))+
  scale_color_viridis_c(trans= "log10", option= "mako", begin = 0.1, end= 0.8, direction= -1,
                        labels= scales::number)+
  theme_classic()+
  labs(x=expression(CO[2]~departure~(mu*mol~L^-1)), y= expression(O[2]~departure~(mu*mol~L^-1)),
       color= expression(Discharge~(m^3~s^-1)))+
  theme(legend.position= "inside", legend.position.inside =  c(.8,.88), legend.key.width = unit(5, "mm"))

#check the plot
plot_ellipse_rcc




#calculate the metrics  by discharge
ellipse_metrics <- sim_out %>%
  filter(date > as.Date("2020-07-13")) %>% #take only one day of the simulation
  mutate(CO2dep_mmolm3 = co2.mod- co2.air,
         O2dep_mmolm3 = o2.mod-o2.air) %>% 
  summarise(meanCO2= mean(CO2dep_mmolm3),
            meanO2 = mean(O2dep_mmolm3), 
            n = n(),
            offset = meanCO2 + meanO2,
            elipse = elipse(CO2dep_mmolm3, O2dep_mmolm3),
            reg = regression_metrics(CO2dep_mmolm3, O2dep_mmolm3), .by= "q.m3s"  ) %>% 
  unnest(c(elipse, reg))


# make the simple plots for the metrics
offset_plot <- ellipse_metrics %>% 
  ggplot(aes(q.m3s, offset ))+
  geom_line(linewidth= 1.5)+
  scale_x_log10()+
  scale_y_continuous(limits = c(0,50), breaks = c(0,25,50))+
  labs(y= "Offset", x = "")+
  theme_classic()

length_plot <- ellipse_metrics %>% 
  ggplot(aes(q.m3s, length ))+
  geom_line(linewidth= 1.5)+
  scale_x_log10()+
  scale_y_continuous(limits = c(0,100), breaks= c(0,50,100))+
  labs(y= "Length", x = "")+
  theme_classic()

eq_plot <- ellipse_metrics %>% 
  ggplot(aes(q.m3s, EQ ))+
  geom_line(linewidth= 1.5)+
  scale_x_log10()+
  scale_y_continuous(limits = c(0,1), breaks = c(0, .5, 1))+
  labs(y= "EQ", x= expression(Q~(m^3~s^-1)))+
  theme_classic()

width_plot <- ellipse_metrics %>% 
  ggplot(aes(q.m3s, width ))+
  geom_line(linewidth= 1.5)+
  scale_x_log10()+
  scale_y_continuous(limits = c(0,14), breaks = c(0,5,10))+
  labs(y= "Width", x= "")+
  theme_classic()

plot_metrics <- length_plot + width_plot + offset_plot + eq_plot  + 
  plot_layout(ncol= 1) 

#put together the final figure
plot_ellipse_rcc + 
  plot_metrics + 
  plot_layout(ncol= 2, widths = c(2.6, 1)) +
  plot_annotation(tag_levels = 'a')

ggsave("plots/main/fig5_ellipses_rcc.png", width = 7, height = 5.5)  




    