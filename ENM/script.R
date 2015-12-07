### EDITADO PARA INCLUIR O SEGUNDO FILTRO, DE PONTOS UNICOS POR PIXEL. 7/12/2015
source("./fct/modelos.R")
source("./fct/final_model.R")
source("./fct/ensemble.R")
library("raster")
library("XML")
library("snowfall")
library("rJava")
predictors <-  stack(list.files("./env",full.names = T)[1])
####Ler os dados de ocorrencia - tabela por registro----
registros<-read.csv("./data/Atlantic_Domain_Ari_RevTaxon_Final.csv",sep=";")

#checando NAs----
names(registros)
#extrai os valores ambientais em caad registro e inclui o numero do pixel
spp.env <- extract(predictors,registros[,c("lon","lat")],cellnumbers=T)
#junta com os registros para ter coordenadas, espécies e ambiente no mesmo
registros2 <- cbind(registros,spp.env)
#cria a coluna de nomes
registros2$sp <- paste(registros2$genero,registros2$especie)

head(registros2)[,1:6]
dim(registros2)

#filtra pelo ambiente, os casos (linhas da tabela ambiental que estejam completos)
spp.env.filt <- registros2[complete.cases(spp.env),]
head(spp.env.filt)
dim(spp.env.filt)

# agora o passo é tirar as linhas que tem espécie igual E ambiente igual
a <- spp.env.filt[,c("sp","cells")]
spp.filt <- spp.env.filt[!duplicated(a),]
head(spp.filt)
dim(spp.filt)

names.sp <- unique(spp.filt$sp)
length(names.sp)
names.sp <- names.sp[table(spp.filt$sp)>10] #considerando apenas spp com mais de 10 registros. AGORA É 10 REGISTROS UNICOS, SEM DUPLICADOS ENTRE PIXEIS
length(names.sp)###O NUMERO DE ESPECIES CON N>10 MUDOU, CAIU PARA 2115 ANTES ERA 2200 E POUCO
occs <- spp.filt[,c("sp","lon","lat")]
dim(occs)

#### ATE AQUI, O RESTO FICOU IGUAL


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
args(dismo.mod)
N <- sample(1:2115,2)
lapply(names.sp[N],
       dismo.mod,
       occs = spp.filt,
       predictors = predictors,
       buffer = T,
       buffer.type ="max",
       maxent = T,
       Bioclim = T,
       Domain = F,
       Mahal = F,
       GLM = T,
       RF = T,
       SVM = T,
       SVM2 = T,
       part = 3,
       seed = 712,
       output.folder = "buffermax")
lapply(names.sp[1:10],final.model,select.partitions=T,weight.partitions = T)
lapply(names.sp[1:10],ensemble)
args(ensemble)

plot(mean(bin.sel))

lapply(names.sp[1:10],dismo.mod,maxent=T,Bioclim=T,Domain=T,Mahal=T,GLM=T,RF=T,SVM=T,
       SVM2=T, part=3, seed=2311,predictors=predictors)
q()