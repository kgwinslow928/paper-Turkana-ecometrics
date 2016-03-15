# Computational experiments for Turkana ecometrics manuscript

This repository provides source code for reproducing the computational analysis reported in the manuscript **An ecometric analysis of the fossil mammal record of the Turkana Basin** by Fortelius et al 2016.

Scripts are numbered in the order of appearance in the manuscript text. 

	run1_fit_modern_regression.R
	
Fits regression models on the modern data, plots Figure 1 (visualization of models) and Figure 10 (modern data maps). The data used for modeling is in.

	data/data_Africa_modern.csv
	
It is a grid cell level data, containinng MAP and MAT estimates, mean HYP and mean LOP for each grid cell. We are not able to store a copy of IUCN occurence data here, the grid-level aggregate data is enough to reproduce the models by running the script above.


	run2_estimates_for_fossils.R
	
Adds precipitation and temperature estimates to fossil datasets. There are several complilations of the fossil datasets required for different parts of the analysis. 

	data/Turkana_47_for_code.csv
	data/Turkana_47_for_code_bins.csv
	data/Turkana_47_for_code_members.csv
	data/Turkana_47_for_code_carnivors.csv
	
The first one is the main dataset at specimen level, including ComLocs that have three identified species or more. 'bins' data is aggregated (unique species identified) at time bin level. 'members' data is aggregated at member level. 'carnivors' data is processed the same way as the main dataset at ComLoc level, but it includes carnivor species, this is needed for reproducing Figure 9 -- proportions of species and specimen. 

	run3_regression_analysis.R	
	
Produces analysis figures: Figure 1,3,4,6,7,9. The plots are stored in the folder called 'figures'.	
	
