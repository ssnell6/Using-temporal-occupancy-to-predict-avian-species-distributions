# SDMs

Data and code for project "Using temporal occupancy to predict avian species distributions"

# The SDMs repo includes four folders:
- Data
- Figures
- maxent
- Scripts

# Data 
- The data folder contains all processed data frames used in analyses, including RMSE output, temporal, and spatial cross-validations.

# Figures
- The figures folder contains all visual output in PDF form.

# maxent
- The maxent folder contains the maxent software required to run analyses.

# Scripts
- The scripts folder includes all scripts used in analyses.
- The abiotic script and abiotic script LL were used to process the raw environmental data in 40 km averages.
- The BBS analysis script filtered and cleaned the raw Breeding Bird Survey data for our anlayses.
- The Figure 1 MS script created the maps for the Yellow-throated Viero, figure 1 in the manuscript.
- The GIS script cleaned all spatial data for analsyes, such as creating the set of absences falling within each species range.
- The raster stack spp script cleaned all raw enviromental data for maxent analyses.
- The SDM analysis script creates the final input data frame for analysis, and runs a loop to calculate RMSE for each method and each species. It also creates the RMSE scatter and density plots, as well as the range occupancy plot.
- The SDM anlaysis 5 year script runs all the same analyses as the SDM analysis script, using 5 years of BBS sampling as input rather than 15 years.
- The SDM analysis maxent presence only script runs all main analyses for maxent, which required slightly different input format in order to run the maxent program.
- The SDM analysis 5 year maxent presence only script runs all main analyses for maxent using 5 years of BBS data rather than 15, which required slightly different input format in order to run the maxent program.
- The SDM crossval script runs all cross-valiation analyses and creates the cross-validation figure in the manuscript.
- The SDM crossval maxent presence only script runs all cross-validation analyses for maxent, which required slightly different input format in order to run the maxent program.







