#MidgeGrowth2020

This contains the data and code associated with "The balance between growth and future resource supply for an aquatic insect feeding on contemporaneously growing algae." 

Completely unaltered experimental data are in the raw data folder and analysis ready code are in the clean data folder. The code used to prepare the datasheets used in the experiment is in scripts/PrepareData.R. This includes code to convert observations of dissolved oxygen to GPP and to assign instars to midge larvae based on their head capsule width. The products of this cleaning are in clean data, with one additional dataset (routine midge head capsules from 2013-2020). These data come from Wetzel et al. (2021) https://doi.org/10.6084/m9.figshare.14095697.v3.

scripts/MG_Functions.R includes the functions we used multiple times and across different scripts. It also included universal aesthetics. This includes the function we optimize over to fit the model, the code used to project the data, the code used to estimate midge weights and to estimate growth rates and most unit conversions. 

scripts/Experiment_Analyses.R includes the statistical analyses of experimental data to address the questions: what effect did initial algal abundance and midge presence have on GPP, midge abundance, and midge size. It includes all analyses for tables 1-3 and the code to generate figure 1.

scripts/Experiment_ProductionAnalyses.R includes the code to assess the relationship between instantaneous measures of GPP and midge growth. This includes the non-parametric bootstrapping procedure, fitting the measurement error model, and producing figure 2.

scripts/MG_Model.R includes the code used to fit the model to the data and to project data under altered maximum resource growth rates and consumption rates. It also includes the code to generate figures 3 and 4.

scripts/Supplement.R includes all the code used in the supplement of the manuscript. 