# R-model backend


This repository holds all the data and functions necessary to perform environmental niche models (ENM) using `dismo` and associated packages both in console and in the graphic interface Model-R.

<!--acho que isto deveria sair, junto com os arquivos - The **ARF_sp_DataPreparation** file has the scripts for data preparation using GDAL tools and Linux terminal. Author: Felipe Barros (f.barros@iis-rio.org), this script was developed to ensure the documentation of all procedures and that all rasters have the same propreties. -->

The folder **ENM** has the following structure:

1. At the root of ENM are the scripts where all the modeling procedures are executed and documented 

<!-- eu manteria apenas script.R -->

2. `/data` holds the species occurrences file, `FLORA_occs_final.csv`. The species occurrences com from the [NeoTropTree](http://prof.icb.ufmg.br/treeatlan/treeatlanE_2_historia.htm) project for plants. 
This is a final file that was already cleaned to adequate the original species names to the current correct names of the Flora of Brazil (16/02/2016). The file columns are the scientific name, latitude and longitude. 

<!-- and the shapefile for the BAF.--> 
<!-- Subfolder `tax flora` documents the steps performed during the scientific name cleaning, to adequate the original species names to the current correct names of the Flora of Brazil (16/02/2016).-->

2. `/env` has the predictor variables. These come from the 19 bioclimatic variables in Worldclim and the six topography variables in CGIAR. They have been transformed into a set of six "eigenvariables", a PCA was performed for the entire area and the first six components were spatialized, retaining 95% of the variation. The resolution of these predictors is 1km x 1km and the final extent is the BAF.

<!-- precisamos trocar essa resolução>
<!-- The function for performing this transformation is in `fct/eigenvariables.fct.R`.-->

3. `/fct` has all the functions needed to perform the ENMs.
    
    + 'modelos.R' defines function`dismo.mod`, that allows to choose which predictors (as a raster::stack object) and which algoritms are used. The algorithm options are maxent, bioclim, domain, mahalanobis distance, glm, random forests and two svm implementations. It also allows to choose the number of partitions (default=3), and a seed for reproducibility. The output folder is where models will be stored. The projection to other datasets is partially implemented. The number of pseudoabsences can also be modified (default=500). `dismo.mod` produces models for every partition but does not join them. 
    
    +   Joining partitions is performed by the function `final_model` in `final_model.R`. Two ways of joining partitions are used: 

        1. Models are selected according to the value of their maximum TSS value, then they are cut but the threshold that maximizes TSS. The final model is a mean of these binary models, cut again to match a majority consensus (at least 50% of the selected models predict each final area). This option is called Final.mean.bin7 due to historical reasons, 
        2. Models are selected according to their maximum TSS value, then a mean is performed and this mean is cut by the mean of the thresholds that maximize each selected model. This option is called Final.bin.mean3.
    
+ The `ensemble` (`ensemble.R`) function receives the species name as parameter, then the general output folder where all model were saved and the output where final models were stored. `Which.models` allows the user to select which type of models are to be averaged (options are Final.bin.mean3 and Final.mean.bin7). The function selects those models and averages them, saving them in the specified `output.folder= (default= "ensemble")`.

Any post-processing operation may be performed manually, such as generating potential richness or checking the variation between models (by calculating the standard deviation between algorithms, for example). 

More information about the file structure and operations performed can be found in the repository [wiki](https://github.com/Model-R/Back-end/wiki). 


<!--precisamos voltar ao final.model original-->



 
