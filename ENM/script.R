### EDITADO PARA INCLUIR O SEGUNDO FILTRO, DE PONTOS UNICOS POR PIXEL. 7/12/2015
### tirei toda a parte de limpeza de dados para outro script

source("./fct/modelos.R")
source("./fct/final_model.R")
source("./fct/ensemble.R")
library("raster")
library("XML")
library("snowfall")
library("rJava")

####eigenvariables
#source("./fct/eigenvariables.fct.R")
#library(raster)
#var<-stack(list.files(path='../1km/', pattern='.asc$', full.names = TRUE))
#eigenvariables.fct(var, '1Km_', .95)
#####


occs <- read.csv("./data/registroslimpos.csv",row.names=1,col.names=c("sp","lon","lat"))
dim(occs)
head(occs)
is(occs)
class(occs$sp)
names.sp <- unique(occs$sp)
names.sp == sort(names.sp)
names.sp <- sort(names.sp)####ISTO AQUI É ESSENCIAL PARA FAZER O TABLE
names(table(occs$sp)) == names.sp
length(names.sp)
sum(table(occs$sp)>10)#todas têm mais de 10 registros


#iniciar snowfall
#     sfInit(parallel=T,cpus=24)
# 	#exporta variáveis e funções e pacotes
# 	sfExportAll()
# 	sfLibrary(rJava) #removido aspas.
# 	sfLibrary(raster)
# 	#sfSource("./fct/modelos.R")
# 	#sfSource("./fct/final_model.R")
# 	sfSource("./fct/ensemble.R")
# 	#executa a função
# 	#tInicial <- Sys.time()
# 	#sfClusterApplyLB(names.sp,dismo.mod, predictors = predictors, Bioclim = T, Domain=T, Mahal=T, GLM=T, RF=T, SVM=T, SVM2=T, maxent=T, occs=spp.filt,seed=123, n.back=1000)
# 	#tFinal <- Sys.time()
# 	#tFinal - tInicial
#
# 	#tInicial <- Sys.time()
# 	#sfClusterApplyLB(names.sp, final.model)
# 	#tFinal <- Sys.time()
# 	#tFinal - tInicial
# 	sfClusterApplyLB(names.sp, ensemble, input.folder1 = "models", input.folder2 = "presfinal")
#
# sfStop()
#write.csv(occs,"registroslimpos.csv")

args(dismo.mod)
set.seed(712)
N <- sample(1:2114,30)
names.sp[N]

lapply(names.sp[N],
       dismo.mod,
       occs = occs,
       predictors = predictors,
       buffer = T,
       buffer.type ="max",
       maxent = T,
       Bioclim = F,
       Domain = F,
       Mahal = F,
       GLM = T,
       RF = T,
       SVM = T,
       SVM2 = F,
       part = 3,
       seed = 712,
       output.folder = "max")

args(final.model)
sort(names.sp[N])

lapply(sort(names.sp[N]),
       final.model,
       input.folder="nobuffer",
       output.folder="presfinal")

args(ensemble)
lapply(sort(names.sp[N])[2:12],
       ensemble,
       input.folder1="nobuffer",
       input.folder2="presfinal",
       occs=occs,
       which.models=c("Final.bin.mean3", "Final.mean.bin7"),
       output.folder="ensemble")


#### teste final model com anfibios

anf.names <- list.files("./mods")
args(finalModel)


lapply(anf.names[2],
       finalModel,
       input.folder="mods",
       output.folder="finalmods")

source("fct/ensemble.R")
args(ensemble)
lapply(anf.names[2],
       ensemble,
       input.folder1="mods",
       input.folder2="finalmods",
       output.folder = "ensemble",
       occs=occs
       )
