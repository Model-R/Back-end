source("./fct/modelos.R")
#source("./fct/final_model.R")
source("./fct/ensemble.R")
library("raster")
library("XML")
library("snowfall")
library("rJava")
predictors <-  stack(list.files("./env",full.names = T)[1])
####Ler os dados de ocorrencia - tabela por registro----
registros<-read.csv("./data/Atlantic_Domain_Ari_RevTaxon_Final.csv",sep=";")
#checando NAs----
spp.filt <- extract(predictors,registros[,6:5])
spp.filt <- registros[complete.cases(spp.filt),]
spp.filt$sp <- paste(spp.filt$genero,spp.filt$especie)
names.sp <- unique(spp.filt$sp)
names.sp <- names.sp[table(spp.filt$sp)>10] #considerando apenas spp com mais de 10 registros
#iniciar snowfall 
sfInit(parallel=T,cpus=24)
	#exporta variáveis e funções e pacotes
	sfExportAll()
	sfLibrary(rJava) #removido aspas.
	sfLibrary(raster)
	#sfSource("./fct/modelos.R")
	#sfSource("./fct/final_model.R")
	sfSource("./fct/ensemble.R")
	#executa a função 
	#tInicial <- Sys.time()
	#sfClusterApplyLB(names.sp,dismo.mod, predictors = predictors, Bioclim = T, Domain=T, Mahal=T, GLM=T, RF=T, SVM=T, SVM2=T, maxent=T, occs=spp.filt,seed=123, n.back=1000)
	#tFinal <- Sys.time()
	#tFinal - tInicial

	#tInicial <- Sys.time()
	#sfClusterApplyLB(names.sp, final.model)
	#tFinal <- Sys.time()
	#tFinal - tInicial
	sfClusterApplyLB(names.sp, ensemble, input.folder1 = "models", input.folder2 = "presfinal")

sfStop()
args(dismo.mod)
dismo.mod(names.sp[1],maxent=T,Bioclim=T,Domain=T,Mahal=T,GLM=T,RF=T,SVM=T,
       SVM2=T, part=3, seed=2311,predictors=predictors)

lapply(names.sp[1:3],dismo.mod,maxent=T,Bioclim=T,Domain=T,Mahal=T,GLM=T,RF=T,SVM=T,
       SVM2=T, part=3, seed=2311,predictors=predictors)
q()