# cca_pls_comparison_ABCD
This repository contains the scripts used in the manuscript titled "Comparing the stability and reproducibility of brain-behaviour relationships found in Canonical Correlation Analysis and Partial Least Squares within the ABCD sample"


The `data_cleaning` folder contains the scripts used to organize and clean the data from the ABCD study to prepare them for use. It contains 3 files which organize the 3 different analysis sets: 

1. The first analysis relating CBCL scores to cortical thickness
2. The second analysis relating NIH scores to cortical thickness
3. The post-hoc analyses linking CBCL scores (among a subset of participants with higher presentation of scores) to cortical thickness

The `Aim1_withinMethod_generalization` folder contains the scripts used to perform the full CCA or PLS analyses. Only one script (for each approach) is included in the folder. However, these scripts have been used for the first, second, and posthoc analyses with the original resiualized matrices changing depending on the analysis being done. 

The `Aim2_betweenMethod_generalization` folder contains 2 folders: 

1. `Functions` 

	The `Functions` folder contains all the resampling analysis functions performed in the current study: permutation, train-test, split-half, bootstrap. This folder contains the scripts used for CCA and PLS separately and only contains the functions that analyzed a Spearman's cross-product matrix given that this was used in the main paper. Each script in the Functions folder produces an Rdata folder containing a list of the resampled coefficients of interest. 

2. `analyzingResults`

	The `analyzingResults` folder contains the scripts which analyze the Rdata file produced by the scripts in the `Functions` folder. 


