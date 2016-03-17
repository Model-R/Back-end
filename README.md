# ARF_spatial_planning


This repository holds all the data and functions necessary to perform environmental niche models (ENM) for species of plants, amphibians and birds from the Brazilian Atlantic Forest (BAF). The resulting ENM will be used in a spatial prioritization analysis to find areas most suitable for ecological restoration in the BAF.

The **ARF_sp_DataPreparation** file has the scripts for data preparation using GDAL tools and Linux terminal. Author: Felipe Barros (f.barros@iis-rio.org), this script was developed to ensure the documentation of all procedures and that all rasters have the same propreties.

The folder **ENM** has the following structure:
 
1. `/data` holds the species occurrences and the shapefile for the BAF. The species occurrences com from the [NeoTropTree](http://prof.icb.ufmg.br/treeatlan/treeatlanE_2_historia.htm) project for plants. Subfolder `tax flora` documents the steps performed during the scientific name cleaning, to adequate the original species names to the current correct names of the Flora of Brazil (16/02/2016).
2. `/env` has the predictor variables. These come from the 19 bioclimatic variables in Worldclim and the six topography variables in CGIAR. They have been transformed into a set of six "eigenvariables", a PCA was performed for the entire area and the first six components were spatialized, retaining 95% of the variation. The function for performing this transformation is in `fct/eigenvariables.fct.R`. The resolution of these predictors is 1km x 1km and the final extent is the BAF.
3. `/fct` has all the functions needed to perform the ENMs.  
    3.1 Function`dismo.mod` allows to choose which predictors (stack) and which algoritms are used. The options are maxent, bioclim, domain, mahalanobis distance, glm, random forests and two svm implementations. It also allows to choose the number of partitions (default=3), and a seed for reproducibility. The output folder is where models will be stored. The projection to other datasets is partially implemented. The number of pseudoabsences can also be modified (default=500).

`dismo.mod` produces models for every partition but does not join them. Joining partitions is performed by the function `final_model`. In `final_model`, two ways of joining partitions are used: 

1. Models are selected according to the value of their maximum TSS value, then they are cut but the threshold that maximizes TSS. The final model is a mean of these binary models, cut again to match a majority consensus (at least 50% of the selected models predict each final area). This option is called Final.mean.bin7 due to historical reasons,
2. Models are selected according to their maximum TSS value, then a mean is performed and this mean is cut by the mean of the thresholds that maximize each selected model. This option is called Final.bin.mean3,

The `ensemble` function receives the species name as parameter, then the general output folder where all model were saved and the output where final models were stored. Which.models allows the user to select which type of models are to be averaged (options are Final.bin.mean3 and Final.mean.bin7). The function selects those models and averages them, saving them in the specified output.folder= (default= "ensemble").

Any further operation may be performed manually, such as generating potential richness or checking the variation between models (by calculating the standard deviation between algorithms, for example). 

Other supplementary functions are also stored in `fct`.

Finally, at the root of ENM are the scripts where all the modeling procedures are executed and documented,




 
