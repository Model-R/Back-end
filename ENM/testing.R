library("raster")
library("XML")
library("rJava")
library("rgdal")
library("rgeos")
library("dismo")

variaveis_preditoras <- stack(list.files("/home/guilherme/repos/ARF_spatial_planning/ENM/env",pattern="gri",full.names = T)[1])

mascara <- readOGR(dsn="/home/guilherme/repos/ARF_spatial_planning/ENM/data", layer="Bioma_MA1148")

coordenadas <- read.csv("/home/guilherme/repos/ARF_spatial_planning/ENM/data/FLORA_occs_final.csv")

especie <- "Zygia latifolia (L.) Fawc. & Rendle"

ocorrencias <- coordenadas[coordenadas$sp == especie, c("lon", "lat")]

source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/modelos.R")

dismo.mod(sp = especie,
	  occs = coordenadas,
	  predictors = variaveis_preditoras,
	  buffer = FALSE,
	  maxent = T,
	  Bioclim = T,
	  Domain = T,
	  Mahal = F,
	  GLM = F,
	  RF = F,
	  SVM = F,
	  SVM2 = F,
	  part = 3,
	  seed = 512,
	  models.dir = "/tmp/antigo",
	  project.model = F,
	  projections = NULL,
	  mask = mascara,
	  n.back = 500)
