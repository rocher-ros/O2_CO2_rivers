## Load packages and files ----
# Load libraries
library(patchwork)
library(scico)
library(ggridges)

#Load all custom functions for the model
source("scripts/model_o2_co2/model_dic_o2.R")

#load metab data from Appling, keep only good days
daily_metab <- read_delim("prepared data/river data/Appling2019/daily_predictions.tsv") %>% 
  filter(GPP>=0, ER <= 0, GPP.Rhat < 1.1, ER.Rhat < 1.1, K600.Rhat < 1.1)

#get site average data for all variables, calculated in script # 4
site_metab <- read_csv("prepared data/river data/Appling2019/site_avgs_gwinputs.tsv") %>% 
  select(site_name, ends_with("_median"), gw_frac, spQ_mmday, tot_area) %>% 
  rename_with(~ str_remove(., "_median"), everything()) %>% 
  mutate(gw_frac = ifelse(gw_frac > 100 | gw_frac < 0, NA, gw_frac))

gw_usgs <- read_csv("prepared data/groundwater data/usgs_gw_co2_o2_site_avg.csv") %>% 
  mutate(carbon_dioxide_mg_l = ifelse(carbon_dioxide_mg_l > 150, NA, carbon_dioxide_mg_l),#extreme pco2 values, romving above 150.000 ppm
         co2_umol_l = carbon_dioxide_mg_l/44*1e+3,
         co2_ppm = co2_umol_l*1e+3* 0.04477565,
         oxygen_mg_l = ifelse(oxygen_mg_l > 13, NA, oxygen_mg_l),
         oxygen_umol_l= oxygen_mg_l/32*1e+3)

#read the water chemistry file, to get alkalinity
chemistry_usgs <- read_csv("prepared data/river data/USGS_data/chemistry_site_avg.csv")

#join to the previous file
met_with_chem <- site_metab %>% 
  mutate(site = str_remove_all(site_name, "nwis_") %>% 
           paste0("USGS-", .)) %>% 
  left_join(chemistry_usgs, by = "site") %>% 
  mutate(alk_mmol_m3 = alkalinity_mg_l_ca_co3*0.02*1000)

#draw some distributions for PQ and RQ
quotients_distr <- tibble(pq = rnorm(10000, mean= 1.18, sd= .31), #from literature values (SM) #1.18, sd= 0.31
                          rq = rnorm(10000, mean = 1.2, sd= 0.45)) # from Berggren et al. 2012


#some stats for the paper
met_with_chem %>% 
  summarise(across(c(GPP, ER, K600), list(min = min, max= max, mean=mean, sd= sd, count = ~ n())))

met_with_chem %>% 
  drop_na(alk_mmol_m3) %>% 
  summarise(across(c(alk_mmol_m3), list(min = min, max= max, mean=mean, sd= sd,  count = ~ n())))

gw_usgs %>% 
  drop_na(co2_umol_l) %>% 
  summarise(across(c(co2_umol_l), list(min = min, max= max, mean=mean, sd= sd, count = ~ n())))

gw_usgs %>% 
  drop_na(oxygen_umol_l) %>% 
  summarise(across(c(oxygen_umol_l), list(min = min, max= max, mean=mean, sd= sd, count = ~ n())))

#now eyeballing co2, in  mmol/m3 for the model,
# to get this in ppm, do
#  co2_ppm*kh*1000, which for 15C temp and 1 atm would be: umol_col*1000* 0.04477565
#gw_distr <- tibble(gw_o2 = rexp(10000, rate = .4),
#                   gw_co2 = rlnorm(10000, mean= 0, sd = 0.5)*500) 



#calculate average values to use as fix values
avg_values <- site_metab %>% 
  summarise(across(c(ER, GPP, K600, discharge, depth, gw_frac), \(x) mean(x, na.rm = TRUE))) %>% 
  mutate(PQ = mean(quotients_distr$pq),
         RQ = mean(quotients_distr$rq),
         alk_mmol_m3 = mean(met_with_chem$alk_mmol_m3, na.rm = TRUE),
         gw.O2conc = mean(gw_usgs$oxygen_mg_l, na.rm = TRUE), #[gw_usgs$monitoring_location_type_name == "Subsurface"]
         gw.CO2conc = mean(gw_usgs$co2_umol_l, na.rm = TRUE))

## Common plotting settings ----


#quantiles to extract data from the each variables
q_selected <- c(0.1,0.3,0.5,0.7,0.9)

#ggplot commands for the ellipse plots
ellipse_plotting_theme <- list( geom_point(size=1),
                                geom_polygon(alpha=.2, linewidth = 0),
                                geom_abline(slope=-1, intercept = 0, linetype= 2),
                                geom_hline(yintercept = 0),
                                geom_vline(xintercept = 0),
                                scale_x_continuous(limits=c(-10, 200)),
                                scale_y_continuous(limits=c(-110, 50)),
                                theme_classic(),
                                theme(legend.position = "none", panel.border =element_rect(linewidth = 1, fill= NA),
                                      plot.title = element_text(hjust=0.5))
                                )

#ggplot commands for the density plots
density_plotting_theme <- list(geom_density_ridges_gradient(linewidth= 0),
                               theme_classic(),
                               theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                                     axis.line.y = element_blank(), axis.text.x= element_text(color= "black"),
                                     plot.background  = element_rect(linewidth = 0.3, fill= "white", color="gray80"),
                                     plot.margin = margin(3,10,-10,-7))
                               )


## ER plot ----

#calculate selected quantiles
ER_qs <- site_metab %>% 
  reframe(ER =quantile(ER, q_selected, na.rm = TRUE),
            q= q_selected)

#run the model with all fixed except ER
pars_ER <- 
  tibble(ER.day = ER_qs$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         bigK = avg_values$K600,
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         PQ = avg_values$PQ,
         RQ = avg_values$RQ,
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_ER <- pars_ER %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_er <- sim_out_ER %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= ER.day*32/1e+3, fill= ER.day*32/1e+3, group= ER.day*32/1e+3))+
  ellipse_plotting_theme +
  scale_color_scico( palette= "lajolla", begin = 0.2, end=0.8)+
  scale_fill_scico( palette= "lajolla", begin = 0.2, end=0.8)+
  labs(x="", y= "", title= expression(ER))

density_er <- 
  ggplot(site_metab, aes(x = ER, y = 1, fill = after_stat(x))) + 
  density_plotting_theme +
  geom_segment(data= ER_qs, aes(x=ER, xend= ER, y=1, yend= 1.4),  color= "gray30", linetype=1)+
  geom_text(data=ER_qs, aes(x=ER, y= 1.6, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico(name= expression(g~O[2]~m^3~L^-1), palette= "lajolla", begin = 0.1, end=0.8)+
  scale_color_scico(name= expression(g~O[2]~m^3~L^-1), palette= "lajolla", begin = 0.1, end=0.8)+
  scale_y_continuous( expand = c(0,0), limits= c(1,1.7))+
  scale_x_continuous(limits=c(quantile(site_metab$ER, 0.05),
                              quantile(site_metab$ER, 0.95)), expand = c(0,0),
                     breaks= c(-2,-4,-6,-8,-10), name= expression(g~O[2]~m^2~d^-1))+
  labs(x= "", y="")+
  guides(x = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x = element_text( hjust = 2, vjust = 7.1,   color= "gray20", size= 10),
        plot.margin = margin(3,10,-12,-7))


plot_er <- ellipses_er + inset_element(density_er, left = 0.01, bottom = 0.01, right = .68, top= .35)


## GPP plot ----
GPP_qs <- site_metab %>% 
  reframe(GPP =quantile(GPP, q_selected),
          q= q_selected)

#run the model with all fixed except GPP
pars_GPP <- tibble(GPP.day = GPP_qs$GPP/32*1e+3, 
                   ER.day = avg_values$ER/32*1e+3, 
                   bigK = avg_values$K600,
                   q.m3s = avg_values$discharge,
                   depth.m= avg_values$depth,
                   alk.mmolm3 = avg_values$alk_mmol_m3, 
                   PQ = avg_values$PQ,
                   RQ = avg_values$RQ, 
                   gw.frac = avg_values$gw_frac/100, 
                   gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
                   gw.CO2conc = avg_values$gw.CO2conc) 



## run the simulation with the set parameters
sim_out_GPP <- pars_GPP %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


#plot of ellipses
ellipses_gpp <- sim_out_GPP %>% 
  filter(date== as.Date("2020-07-07")) %>% 
ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= GPP.day*32/1e+3, fill= GPP.day*32/1e+3, group= GPP.day*32/1e+3))+
 # annotate("text", x= 174, y= -105, label= expression(g~O[2]~m^2~d^-1), parse= TRUE)+
  ellipse_plotting_theme +
  scale_color_scico( palette= "bamako",direction=-1, begin = 0.3, end=0.9)+
  scale_fill_scico(palette= "bamako",direction=-1, begin = 0.3, end=0.9)+
  labs(x="", y= expression(O[2]~departure~(mu*mol~L^-1)), title= expression(GPP))


#density plot as a legend for each variable
 density_gpp <- ggplot(site_metab, aes(x = GPP, y = 1, fill = after_stat(x))) + 
   density_plotting_theme +
  geom_segment(data= GPP_qs, aes(x=GPP, xend= GPP, y=1, yend= 1.7),  color= "gray30", linetype=1)+
  geom_text(data=GPP_qs, aes(x=GPP, y= 1.9, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico( palette= "bamako",direction=-1, begin = 0.2, end=0.9)+
  scale_x_continuous(limits=c(quantile(site_metab$GPP, 0.05),
                              quantile(site_metab$GPP, 0.95)), expand = c(0,0),
                     name= expression(g~O[2]~m^2~d^-1))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,2.1))+
   guides(x = guide_axis(minor.ticks = TRUE))+
   theme(axis.title.x = element_text( hjust = 2, vjust = 7.1,   color= "gray20", size= 10),
         plot.margin = margin(3,10,-12,-7))
  
# put both together
 plot_gpp <- ellipses_gpp + inset_element(density_gpp, left = 0.01, bottom = 0.01, right = .68, top= .35)

 
## K plot ----
K_qs <- site_metab %>% 
  reframe(K600 =quantile(K600, q_selected),
          q= q_selected)

#run the model with all fixed except K
pars_K <- 
  tibble(bigK = K_qs$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         PQ = avg_values$PQ,
         RQ = avg_values$RQ, 
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_K <- pars_K %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_K <- sim_out_K %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x = co2.mod-co2.air, y = o2.mod-o2.air, color = K600.d, fill= K600.d, group= K600.d))+
  ellipse_plotting_theme +
  scale_color_scico(palette= "oslo", begin = 0.3, end = 0.9)+
  scale_fill_scico(palette= "oslo", begin = 0.3, end = 0.9)+
  labs(x="", y= "", title= expr(K[600]))



density_k <- ggplot(site_metab, aes(x = K600, y = 1, fill = after_stat(x))) + 
  density_plotting_theme +
  geom_segment(data= K_qs, aes(x=K600, xend= K600, y=1, yend= 1.2),  color= "gray30", linetype=1)+
  geom_text(data=K_qs, aes(x=K600, y= 1.25, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico(palette= "oslo", begin = 0.3, end=0.97)+
  scale_x_continuous(limits=c(quantile(site_metab$K600, 0.05),
                              quantile(site_metab$K600, 0.95)), expand = c(0,0),
                     name= expression(d^-1))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,1.3))+
  guides(x = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x = element_text( hjust = 1.13, vjust = 7.1,   color= "gray20", size= 10),
        plot.margin = margin(3,10,-12,-7))

plot_k <- ellipses_K + inset_element(density_k, left = 0.01, bottom = 0.01, right = .68, top= .35)


## Plot Alk ----

alk_qs <- met_with_chem %>% 
  reframe(alk_mmol_m3 =quantile(alk_mmol_m3, q_selected, na.rm = TRUE),
          q= q_selected)

#run the model with all fixed except K
pars_alk <- 
  tibble(bigK = avg_values$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         PQ = avg_values$PQ,
         RQ = avg_values$RQ, 
         alk.mmolm3 = alk_qs$alk_mmol_m3, 
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_alk <- pars_alk %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_alk <- sim_out_alk %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= alk.mmolm3, fill= alk.mmolm3, group= alk.mmolm3))+
  ellipse_plotting_theme +
  scale_color_scico(palette= "davos", begin = 0.2, end=0.8)+
  scale_fill_scico(palette= "davos", begin = 0.2, end=0.8)+
  labs(x="", 
       y= "", title= expression("Alkalinity"))

density_alk <- met_with_chem %>% 
  filter(alk_mmol_m3> 0 ) %>% 
  ggplot( aes(x = alk_mmol_m3, y = 1, fill = after_stat(x))) + 
  density_plotting_theme+
  geom_segment(data= alk_qs, aes(x=alk_mmol_m3, xend= alk_mmol_m3, y=1, yend= 1.00057), color= "gray30", linetype=1)+
  geom_text(data=alk_qs, aes(x=alk_mmol_m3, y= 1.00072, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico( palette= "davos", begin = 0.1, end=0.9)+
  scale_x_continuous(limits=c(quantile(met_with_chem$alk_mmol_m3, 0.05, na.rm= TRUE),
                              quantile(met_with_chem$alk_mmol_m3, 0.95, na.rm= TRUE)), expand = c(0,0),
                     breaks= c(500,2500,4500), name= expression(mu*mol~L^-1), minor_breaks = scales::breaks_width(500) )+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,1.0009))+
  guides(x = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x = element_text( hjust = 1.65, vjust = 7.2,  color= "gray20", size= 10),
        plot.margin = margin(3,10,-15,-7))

plot_alk <- ellipses_alk + inset_element(density_alk, left = 0.01, bottom = 0.01, right = .68, top= .35)

## Plot PQ ----

pq_qs <- quotients_distr %>% 
  reframe(pq =quantile(pq, c(0.1,0.3,0.5,0.7,0.9), na.rm = TRUE),
          q= q_selected)

#run the model with all fixed except K
pars_pq <- 
  tibble(bigK = avg_values$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         PQ = pq_qs$pq,
         RQ = avg_values$RQ, 
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_pq <- pars_pq %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_pq <- sim_out_pq %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= PQ, fill= PQ, group= PQ))+
  ellipse_plotting_theme +
  scale_color_scico(palette= "navia", begin = 0.2, end=0.8)+
  scale_fill_scico( palette= "navia", begin = 0.2, end=0.8)+
  labs(x="", y= expression(O[2]~departure~(mu*mol~L^-1)), title= "PQ")


density_pq <- 
  ggplot(quotients_distr, aes(x = pq, y = 1, fill = after_stat(x))) + 
  density_plotting_theme+
  geom_segment(data= pq_qs, aes(x=pq, xend= pq, y=1, yend= 5.3),  color= "gray30", linetype=1)+
  geom_text(data=pq_qs, aes(x=pq, y= 6.3, label=q ), angle = 90, size= 2.6, color= "gray10")+
  scale_fill_scico( palette= "navia", begin = 0.1, end=0.9)+
  scale_x_continuous(limits=c(quantile(quotients_distr$pq, 0.05, na.rm= TRUE),
                              quantile(quotients_distr$pq, 0.95, na.rm= TRUE)), 
                     name="", expand = c(0,0), breaks= c(1,1.2, 1.4))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,7.4))+
  guides(x = guide_axis(minor.ticks = TRUE))

plot_pq <- ellipses_pq + inset_element(density_pq, left = 0.01, bottom = 0.01, right = .68, top= .35)


## Plot RQ ----


rq_qs <- quotients_distr %>% 
  reframe(rq =quantile(rq, q_selected, na.rm = TRUE),
          q= q_selected)

#run the model with all fixed except K
pars_rq <- 
  tibble(bigK = avg_values$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         PQ = avg_values$PQ,
         RQ = rq_qs$rq,
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_rq <- pars_rq %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)

ellipses_rq <- sim_out_rq %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= RQ, fill= RQ, group= RQ))+
  ellipse_plotting_theme +
  scale_color_scico( palette= "bilbao", begin = 0.2, end=0.8)+
  scale_fill_scico( palette= "bilbao", begin = 0.2, end=0.8)+
  labs(y= "", x= "", title= "RQ")

density_rq <- 
  ggplot(quotients_distr, aes(x = rq, y = 1, fill = after_stat(x))) + 
  density_plotting_theme+
  geom_segment(data= rq_qs, aes(x=rq, xend= rq, y=1, yend= 2.9), color= "gray30", linetype=1)+
  geom_text(data=rq_qs, aes(x=rq, y= 3.4, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico(palette= "bilbao", begin = 0.1, end=0.9)+
  scale_x_continuous(limits=c(quantile(quotients_distr$rq, 0.05, na.rm= TRUE),
                              quantile(quotients_distr$rq, 0.95, na.rm= TRUE)), name="", expand = c(0,0))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,4))+
  guides(x = guide_axis(minor.ticks = TRUE))

plot_rq <- ellipses_rq + inset_element(density_rq, left = 0.01, bottom = 0.01, right = .68, top= .35)

## Plot %GW ----

gw_qs <- site_metab %>% 
  reframe(gw_frac =quantile(gw_frac, q_selected, na.rm = TRUE)/100,
          q= q_selected)

#run the model with all fixed except K
pars_gw <- 
  tibble(bigK = avg_values$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         PQ = avg_values$PQ,
         RQ = avg_values$RQ,
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         gw.frac = gw_qs$gw_frac, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_gw <- pars_gw %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_gw <- 
  sim_out_gw %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= log(gw.frac*100), fill= log(gw.frac*100), group= log(gw.frac*100)))+
  #annotate("text", x= 150, y= -105, label= "%")+
  ellipse_plotting_theme +
  scale_color_scico( palette= "glasgow", direction = -1, begin = 0.2, end=0.9)+
  scale_fill_scico( palette= "glasgow",  direction = -1, begin = 0.2, end=0.9)+
  labs(x= expression(CO[2]~departure~(mu*mol~L^-1)), y= expression(O[2]~departure~(mu*mol~L^-1)), title= "Groundwater discharge")


density_gw <- 
  ggplot(site_metab, aes(x = gw_frac, y = 1, fill = after_stat(x))) + 
  density_plotting_theme+
  geom_segment(data= gw_qs, aes(x=gw_frac*100, xend= gw_frac*100, y=1, yend= 1.9), color= "gray30", linetype=1)+
  geom_text(data=gw_qs, aes(x=gw_frac*100, y= 2.1, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico(palette= "glasgow", direction = -1, begin = 0.1, end=1, trans= "log10", name= "")+
  scale_x_log10(limits=c(quantile(site_metab$gw_frac, 0.05, na.rm= TRUE),
                              quantile(site_metab$gw_frac, 0.95, na.rm= TRUE)), name="%", expand = c(0,0),
                guide = guide_axis_logticks(long = 2, mid = 1, short = 0.5))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,2.4))+
  theme(axis.title.x = element_text( hjust = 1.05, vjust = 5.7,  color= "gray20", size= 10),
        plot.margin = margin(3,10,-10,-7))

plot_gw <- ellipses_gw + inset_element(density_gw, left = 0.01, bottom = 0.01, right = .68, top= .35)



## Plot gw CO2 conc----

gw_co2_qs <- gw_usgs %>% 
  reframe(gw_co2 =quantile(co2_umol_l, q_selected, na.rm = TRUE),
          q= q_selected)

#run the model with all fixed except K
pars_gw_co2 <- 
  tibble(bigK = avg_values$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         PQ = avg_values$PQ,
         RQ = avg_values$RQ,
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = avg_values$gw.O2conc/32*1e+3, 
         gw.CO2conc = gw_co2_qs$gw_co2) 

## run the simulation with the set parameters
sim_out_gw_co2 <- pars_gw_co2 %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_gw_co2 <- 
  sim_out_gw_co2 %>% 
  mutate(gw_pco2= (gw.DICconc-alk.mmolm3)) %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= gw_pco2, fill= gw_pco2, group= gw_pco2))+
  ellipse_plotting_theme +
  scale_color_scico( palette= "tokyo", direction = -1, begin = 0.2, end=0.9)+
  scale_fill_scico( palette= "tokyo",  direction = -1, begin = 0.2, end=0.9)+
  labs(y= "", x= expression(CO[2]~departure~(mu*mol~L^-1)), title= expression(CO[2]~"in groundwater"))


density_gw_co2 <- 
  ggplot(gw_usgs, aes(x = co2_umol_l, y = 1, fill = after_stat(x))) + 
  density_plotting_theme+
  geom_segment(data= gw_co2_qs, aes(x=gw_co2, xend= gw_co2, y=1, yend= 1.0016), color= "gray30", linetype=1)+ #from umol/L to ppm *1000*0.04477565
  geom_text(data=gw_co2_qs, aes(x=gw_co2, y= 1.0019, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico(palette= "tokyo", direction = -1, begin = 0.1, end=1)+
  scale_x_continuous(limits=c(quantile(gw_usgs$co2_umol_l, 0.05, na.rm= TRUE),
                         quantile(gw_usgs$co2_umol_l, 0.95, na.rm= TRUE)), breaks = c(250, 1000, 1750),
                     name=expression(mu*mol~L^-1), expand = c(0,0))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,1.0024))+
  guides(x = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x = element_text( hjust = 1.55, vjust = 7.1,  color= "gray20", size= 10),
        plot.margin = margin(3,10,-15,-7))

plot_gw_co2 <- ellipses_gw_co2 + inset_element(density_gw_co2, left = 0.01, bottom = 0.01, right = .68, top= .35)


## Plot gw O2 conc----

gw_o2_qs <- gw_usgs %>% 
  #filter(monitoring_location_type_name == "Subsurface") %>% 
  reframe(gw_o2 =quantile(oxygen_mg_l/32*1e+3, q_selected, na.rm = TRUE),
          q= q_selected)

#run the model with all fixed except K
pars_gw_o2 <- 
  tibble(bigK = avg_values$K600,
         ER.day = avg_values$ER/32*1e+3, 
         GPP.day = avg_values$GPP/32*1e+3, 
         q.m3s = avg_values$discharge,
         depth.m= avg_values$depth,
         PQ = avg_values$PQ,
         RQ = avg_values$RQ,
         alk.mmolm3 = avg_values$alk_mmol_m3, 
         gw.frac = avg_values$gw_frac/100, 
         gw.O2conc = gw_o2_qs$gw_o2, 
         gw.CO2conc = avg_values$gw.CO2conc) 

## run the simulation with the set parameters
sim_out_gw_o2 <- pars_gw_o2 %>% 
  group_nest(row_number()) %>% 
  pull(data) %>% 
  map_dfr(co2_o2_sim)


ellipses_gw_o2 <- 
  sim_out_gw_o2 %>% 
  filter(date== as.Date("2020-07-07")) %>% 
  ggplot(aes(x=co2.mod-co2.air, y=o2.mod- o2.air, color= gw.O2conc, fill= gw.O2conc, group= gw.O2conc))+
  ellipse_plotting_theme +
  scale_color_scico( palette= "lipari", direction = -1, begin = 0.2, end=0.9)+
  scale_fill_scico( palette= "lipari",  direction = -1, begin = 0.2, end=0.9)+
  labs(y= "", x= expression(CO[2]~departure~(mu*mol~L^-1)), title= expression(O[2]~"in groundwater")) 

density_gw_o2 <- 
  ggplot(gw_usgs, aes(x = oxygen_mg_l/32*1e+3, y = 1, fill = after_stat(x))) + 
  density_plotting_theme+
  geom_segment(data= gw_o2_qs, aes(x=gw_o2, xend= gw_o2, y=1, yend= 1.0069), color= "gray30", linetype=1)+
  geom_text(data=gw_o2_qs, aes(x=gw_o2, y= 1.008, label=q ), angle = 90, size= 2.6)+
  scale_fill_scico(palette= "lipari", direction = -1, begin = 0.1, end=1)+
  scale_x_continuous(limits=c(quantile(gw_usgs$oxygen_mg_l/32*1e+3, 0.05, na.rm= TRUE),
                              quantile(gw_usgs$oxygen_mg_l/32*1e+3, 0.95, na.rm= TRUE)), 
                     name=expression(mu*mol~L^-1), expand = c(0,0), #breaks=c(2,6,10),
                     minor_breaks = scales::breaks_width(50))+
  scale_y_continuous( expand = c(0,0), name="", limits= c(1,1.01))+
  guides(x = guide_axis(minor.ticks = TRUE))+
  theme(axis.title.x = element_text( hjust = 1.65, vjust = 7.1,  color= "gray20", size= 10),
        plot.margin = margin(3,10,-15,-7))

plot_gw_o2 <- ellipses_gw_o2 + inset_element(density_gw_o2, left = 0.01, bottom = 0.01, right = .68, top= .35)

## Compose and export all plots ----
plots_all <- 
  plot_gpp + plot_er + plot_k + 
  plot_pq + plot_rq + plot_alk +
  plot_gw + plot_gw_co2 + plot_gw_o2 +
  plot_layout(ncol= 3) 

ggsave("plots/fig2_main_drivers.png", plots_all, height = 8, width = 8, scale=1.22)

ggsave("plots/fig2_main_drivers.pdf", plots_all, height = 8, width = 8, scale=1.2)

ggsave("plots/fig2_main_drivers.svg", plots_all, height = 8, width = 8, scale=1.22)



## Figure to show the time series of physical drivers for the SM ####

dat_example <- sim_out_gw_o2 %>% 
  filter(gw.O2conc < 15.9, date.time > as.POSIXct("2020-07-12 06:00:00")) #get only one level of experiment and three days

plot_light <- 
ggplot(dat_example)+
  geom_ribbon(aes(date.time, ymin= 0, ymax= radiation), fill= "gold")+
  scale_x_datetime(breaks = "12 hours", date_labels = "%H")+
  theme_classic()+
  labs(x="", y= expression(Radiation~(W~m^-2)))

plot_temp <- 
  ggplot()+
  geom_point(data= dat_example, aes(date.time, temp.air), color= "lightblue")+
  geom_point(data= dat_example, aes(date.time, temp.water), color= "blue2")+
  geom_text(aes(x= c(as.POSIXct("2020-07-14 20:00:00"), as.POSIXct("2020-07-14 20:00:00")),
                y= c(15.5, 16.4),
                label= c("Air", "Water")),
            color= c("lightblue", "blue2"), fontface= "bold", hjust= 0)+
  scale_x_datetime(breaks = "12 hours", date_labels = "%H")+
  theme_classic()+
  labs(x="Hour", y= "Temperature (\u00B0C)")

plot_light +
  plot_temp +
  plot_layout(ncol = 1)

ggsave("plots/SM/example_stream.png", scale= .7)

