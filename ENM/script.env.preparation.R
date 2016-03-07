# Script with some changes done in the environmental data
#Reading aspect
library(raster)
aspect1k <- raster(list.files("./1km", full.names=TRUE, pattern='aspect.asc$'))
recs<-which(values(aspect1k==-1))
values(aspect1k)[recs]<-trunc(runif(length(recs),0,36000))
aspect1k<-aspect1k*pi/18000
#plot(aspect1k)

writeRaster(aspect1k, filename = './1km/h_aspect.asc',format='ascii',overwrite=TRUE)
sin.aspect1k<-sin(aspect1k)
writeRaster(sin.aspect1k, filename = './1km/sin.aspect.asc',format='ascii',overwrite=TRUE)
cos.aspect1k<-cos(aspect1k)
writeRaster(cos.aspect1k, filename = './1km/cos.aspect.asc',format='ascii',overwrite=TRUE)

# Renaming variables not be used on PCA analysis:
filesToRename <- (c('./1km/Alt.asc','./1km/h_aspect.asc', './1km/h_flowacc.asc', './1km/h_flowdir.asc'))
class(filesToRename)
for (a in 1:length(filesToRename)){
  file.rename(from = filesToRename[a], to = paste0(filesToRename[a],'XXXX'))
}

# PCA analysis to build eigenvariables
# Issue #1 sourcing from github account I found an Error due its a private account (I imagine) (Felipe)
#library(devtools)
#eigenvariables.fct.R <- "https://github.com/FelipeSBarros-IIS/ARF_spatial_planning/blob/master/ENM/fct/eigenvariables.fct.R"
#source_url(eigenvariables.fct.R)
source("./fct/eigenvariables.fct.R")
library(raster)
var<-stack(list.files(path='../1km/', pattern='.asc$', full.names = TRUE))
eigenvariables.fct(var, '1Km_', .95)