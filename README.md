# Analysis for the article "Emerging patterns of CO2:O2 dynamics in rivers and their link to ecosystem carbon processing"


This repository contains the code to explore how CO2 and O2 can vary in rivers, and what are the relative role of multiple drivers. The article that builds on this analysis is currently in preparation
This project is structured in two parts:
1. A simple model that reproduces CO2 and O2 dynamics in rivers, affected by metabolism, gas exchange and inputs from groudnwater.
2. A data synthesis to constrain and parametrize the model.

The  `scripts/` folder contains three sub-folders:
- `model_o2_co2`, which contains the scripts to create an artificial stream and the model implementation
- `data_cleaning`, which contains the scripts to do some pre-processing of the data
- `article_figures`, containing the code to clean data, run the model, and produce the figures of the article.

The data for the main analysis is deposited in <link-to-repo>, once extracted in a folder called `prepared data`, is possible to run the code.
Below I detail the content on each folder.


## `model_o2_co2`
- script 1...
- 


![Figure 2 of the paper](https://github.com/rocher-ros/O2_CO2_rivers/blob/main/plots/main/fig2_main_drivers.png)
Figure 2 of the paper, showing the relative effect of each model variable on CO2:O2 departure plots. In each panel, the inset shows the kernel distribution function of the empirical data of each variable used in this study (tails of distribution clipped to 0.05-0.95 quantiles of the data, see methods for data sources). For each variable, the model is run at five different quantiles of the data distribution (0.1, 0.3, 0.5, 0.7, 0.9 quantiles, shown in the inset), and the model output as CO2:O2 concentrations is shown in the main plot, with different colours to relate to each quantile used. In each panel, the rest of the variables were kept fixed at their mean values, which were: GPP = 2.56 g O2 m-2 d-1, ER = -4.73 g O2 m-2 d-1, K600 = 7.58 d-1, PQ = 1.19, RQ = 1.19, Alkalinity = 2221 µmol L-1, fraction of groundwater discharge = 1.18 %, CO2 in the groundwater 772 µmol L-1, and O2 in the groundwater:  173 µmol L-1.

