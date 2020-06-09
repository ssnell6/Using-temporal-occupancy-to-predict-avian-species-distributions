# SDMs

Data and code for project "Using temporal occupancy to predict avian species distributions"

The SDMs repository includes four folders and the R project used for all anlayses:
- Data
- Figures
- maxent
- Scripts
- SDMs.Rproj

# Data 
- The data folder contains all processed data frames used in analyses, including RMSE output, temporal, and spatial cross-validations.
- Notable files include:
    - all_expected_pres is the full set of presences and expected absences for all species included in our analysis, created in the GIS script and filtered and cleaned in the SDM analysis script.
    - auc_df is the final RMSE output for each species and method from the SDM analysis script. This file was used to make Figure 2
    - temporal_crossval_df is the output of the temporal cross-validation analyses. The number suffix indicates the threshold used to indicate predicted presences (0.25, 0.5, 0.75).
    - space_cval_rmse is the output of the k-fold spatial cross-validation analysis. The number suffix indicates the threshold used to indicate predicted presences (0.25, 0.5, 0.75).

# Figures
- The figures folder contains all visual output in PDF form resported in the manuscript. Supplemental or old files are in the scratch subfolder.

# maxent
- The maxent folder contains the maxent software required to run analyses.

# Scripts
- The scripts folder includes all scripts used in analyses.
- The abiotic script and abiotic script LL were used to process the raw environmental data in 40 km averages.
- The BBS analysis script filtered and cleaned the raw Breeding Bird Survey data for our anlayses.
- The Figure 1 MS script created the maps for the Yellow-throated Viero, figure 1 in the manuscript.
- The GIS script cleaned all spatial data for analyses, such as creating the set of absences falling within each species range.
- The raster stack spp script cleaned all raw enviromental data for maxent analyses.
- The SDM analysis script creates the final input data frame for analysis, and runs a loop to calculate RMSE for each method and each species. It also creates the RMSE scatter and density plots, as well as the range occupancy plot (Figures 2 and 4).
- The SDM analysis 5 year script runs all the same analyses as the SDM analysis script, using 5 years of BBS sampling as input rather than 15 years, and creates Figure 3.
- The SDM analysis maxent presence only script runs all main analyses for maxent, which required slightly different input format in order to run the maxent program.
- The SDM analysis 5 year maxent presence only script runs all main analyses for maxent using 5 years of BBS data rather than 15, which required slightly different input format in order to run the maxent program.
- The SDM crossval script runs all cross-valiation analyses and creates the cross-validation figure in the manuscript and creates Figure 5 and supplemental figures.
- The SDM crossval maxent presence only script runs all cross-validation analyses for maxent, which required slightly different input format in order to run the maxent program.







