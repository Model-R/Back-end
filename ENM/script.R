# Loading libraries:
library("raster")
library("XML")
library("snowfall")
library("rJava")
library("rgdal")

setwd("./ENM") #porque estou abrindo o projeto inteiro, não a pasta. não sei se seja necessário resolver
# Loading environmental data, study area mask
predictors <-  stack(list.files("./env", full.names = T)[1])
# Cortando pela MataAtlantica:
mascara <- readOGR(dsn = "./data", layer = "Bioma_MA1148")
plot(mascara)
# FLORA ----
# Loading occorrences:
occs <- read.csv("./data/dados_clean.csv")
head(occs)

#só para ter tranquilidade sobre isso, vamos tirar os autores
library(flora)
nomes <- vector()
for (i in seq_along(occs$sp)) {
    nomes[i] <- remove.authors(as.character(occs$sp[i]))
}
head(nomes)
occs$nomes <- nomes
names.sp <- unique(occs$nomes)
length(names.sp)

set.seed(712)

#não usamos mais as funções locais!!!! só para testar que esteja dando igual.
#instalamos o pacote (só é preciso uma vez)
library(devtools)
#devtools::install_github("Model-R/modelr_pkg")
#chamamos o pacote
library(modelr)

# MODELOS.R----
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_144.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
library(rJava)


####Un loop con todas las especies y maxent. Va a sobreescribir lo anterior
lapply(names.sp,
       modelr::do_maxent,
       coordinates = occs[,c(2,3)],
       partitions = 3,
       buffer = "max",
       seed = 310,
       predictors = predictors,
       models.dir = "models",
       n.back = 5000,
       project.model = F,
       mask = mascara)

tif.test <- raster("./models/Abarema langsdorffii/maxent_cont_Abarema langsdorffii_1.tif")
plot(tif.test)

# perfeição <3

lapply(names.sp,
       modelr::do_randomForest,
       coordinates = occs[,c(2,3)],
       partitions = 3,
       buffer = "max",
       seed = 310,
       predictors = predictors,
       models.dir = "models",
       n.back = 5000,
       project.model = F,
       mask = mascara)

lapply(names.sp,
       modelr::do_svm,
       coordinates = occs[,c(2,3)],
       partitions = 3,
       buffer = "max",
       seed = 310,
       predictors = predictors,
       models.dir = "models",
       n.back = 5000,
       project.model = F,
       mask = mascara)

# finalModel ----
args(modelr::finalModel)
lapply(names.sp, modelr::finalModel,
       algoritmos = c("rf", "maxent", "svm"),
       models.dir = "models")

# Ensembles
args(modelr::ensemble)

lapply(names.sp,
       modelr::ensemble,
       occs = occs[,c(2,3)],
       models.dir = "models",
       which.models = c("Final.mean.bin7"),
       consensus = T,
       consensus.level = 0.5)

#Até aqui deve rodar tudo tranquilo. Arquivos mínimos para testar coisas. Eliminei o trecho referente a aves, anfibios e as funções acessórias do Felipe.
