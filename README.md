# Exploring O2 and CO2 dynamics in rivers using a simple ecosystem model
This repository contains the code to explore how CO2 and O2 can vary in rivers, and what are the relative role of multiple drivers. 
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


![Figure 2 of the paper](https://github.com/rocher-ros/O2_CO2_rivers/blob/main/plots/main/fig2_main_drivers.png)
