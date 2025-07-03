# Code of the article "Emerging patterns of CO2:O2 dynamics in rivers and their link to ecosystem carbon processing"

This repository contains the code of the model, as well as to reproduce the results and figures, of the article "Emerging patterns of CO2:O2 dynamics in rivers and their link to ecosystem carbon processing", currently in review in Limnology and Oceanography Letters. This article explores how CO2 and O2 concentrations vary in rivers in response to some of the key ecosystem processes. This is done using some data synthesis to constrain key parameters, as well as a simple model to reproduce CO2 and O2 dynamics based on aerobic metabolism, gas exchange and inputs from groudnwater.

The `scripts/` folder contains three sub-folders: - `model_o2_co2`, which contains the scripts to create an artificial stream and the model implementation - `data_cleaning`, which contains the scripts to do some pre-processing of the data - `article_figures`, containing the code to clean data, run the model, and produce the figures of the article.

The data for the main analysis is deposited in Zenodo (https://doi.org/10.5281/zenodo.13990673). To automatically extract download the data and extract the files in the right folder called `prepared data`, simply run the script `data_setup.R` first. Below I detail the content on each folder.

## `/model_o2_co2`
- `auxiliary_functions.R`: This script has the necessary functions to setup the model.
- `model_dic_o2.R`: The main script with the model. It does not run here, it happens in other scripts.

## `/article_figures`
- `1_lakes_rivers.R`: Script used to collect and compare patterns in CO2:O2 in rivers and lakes, base don the compliation from Vachon et al. 2020 and this study. It produces figure 1.
- `2_model_variables.R`: In this script I use synthesis for different parameters to explore how they afect CO2:O2 dynamics in streams. It produces figure 2.
- `2b_model_variables_mean.R`: This is the same as script 2 but using mean values instead of medians, something that arised during the review process. igure produced is in the SM.
- `3_river_regimes.R`: Script that compares sites with different metabolic regimes, in terms of CO2 and O2 departures. It produces figure 3.
- `4_RiverContiuum.R`: Script that explores patterns of CO2:O2 along the river contiuum. Produces figure 4 and 5.
- `5_SM.R`: This script produces all analysis contained in the supplementary materials of the paper.

## `/data_cleaning`
These scripts need some raw_raw files, available in other publications or the corresponding author if you really want to do everything from scratch.
- `1_stream_prepare_data.R`: Script to fix the stream and rievr data
- `2_final_files.R`: FInal script to prepare all files used in the analysis




![Figure 2 of the paper](https://github.com/rocher-ros/O2_CO2_rivers/blob/main/plots/main/fig2_main_drivers.png) Figure 2 of the paper, showing the relative effect of each model variable on CO2:O2 departure plots. In each panel, the inset shows the kernel distribution function of the empirical data of each variable used in this study (tails of distribution clipped to 0.05-0.95 quantiles of the data, see methods for data sources). For each variable, the model is run at five different quantiles of the data distribution (0.1, 0.3, 0.5, 0.7, 0.9 quantiles, shown in the inset), and the model output as CO2:O2 concentrations is shown in the main plot, with different colours to relate to each quantile used. In each panel, the rest of the variables were kept fixed at their mean values, which were: GPP = 2.56 g O2 m-2 d-1, ER = -4.73 g O2 m-2 d-1, K600 = 7.58 d-1, PQ = 1.19, RQ = 1.19, Alkalinity = 2221 µmol L-1, fraction of groundwater discharge = 1.18 %, CO2 in the groundwater 772 µmol L-1, and O2 in the groundwater: 173 µmol L-1.