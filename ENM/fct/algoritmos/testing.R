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

source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/maxent.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/domain.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/bioclim.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/random_forest.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/glm.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/svm2.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/mahal.R")
source("/home/guilherme/repos/ARF_spatial_planning/ENM/fct/algoritmos/svm.R")


do_domain(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_maxent(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_bioclim(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_randomForest(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_GLM(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_SVM2(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_mahal(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
do_SVM(especie, coordinates = ocorrencias, partitions = 3, buffer = FALSE, seed = 512, predictors = variaveis_preditoras, models.dir = "/tmp/novo", project.model = F, projections = NULL, mask = mascara, n.back = 500)
