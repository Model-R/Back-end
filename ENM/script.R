### EDITADO PARA INCLUIR O SEGUNDO FILTRO, DE PONTOS UNICOS POR PIXEL. 7/12/2015
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

predictors <-  stack(list.files("./env",pattern="1K",full.names = T)[1])
####Ler os dados de ocorrencia - tabela por registro----

###RODEI ISTO PARA CRIAR A NOVA TABELA DE REGISTROS COM OS NOMES CORRIGIDOS. AGORA PODE PULAR PARA A SEGUINTE PARTE

#  registros<-read.csv("./data/Atlantic_Domain_Ari_RevTaxon_Final.csv",sep=";")
#  nomes <- read.delim("./data/lista final.txt")
#
# head(registros)
# registros$sp <- paste(registros$genero,registros$especie)
#
# names(nomes)
# nomes$sp <- paste(nomes$genero,nomes$epiteto)
# length(nomes$sp)
# length(unique(registros$sp))
# setdiff(unique(nomes$sp),unique(registros$sp))
# setdiff(unique(registros$sp),unique(nomes$sp))
#
# registros.final <-  vector(length= length(registros$sp))
# for (i in 1:length(registros.final)) registros.final[i] <- which(registros$sp[i]==nomes$sp,arr.ind=TRUE)
# registros$FINAL <- nomes$FINAL[registros.final]
# unique(registros$FINAL)
# unique(nomes$FINAL)
# names(registros)
#
# library(stringr)
#  spl <- str_split(registros$FINAL,pattern = " ")
#  library(data.table)
#  gen.final <- vector(length=length(registros$FINAL))
#  for (i in 1:length(gen.final)) gen.final[i] <- spl[[i]][1]
# #
#  epithet.final <- vector(length = length(registros$FINAL))
#  for (i in 1:length(gen.final)) epithet.final[i] <- spl[[i]][2]
# #
#  registros$sp.final <- paste(gen.final,epithet.final)
# length(registros$sp.final)
#  rm(epithet.final)
# rm(gen.final)
# rm(spl)
#  #write.csv(registros,file = "./data/registrosnomesfinais.csv")
#
# #registros <- read.csv("./data/registrosnomesfinais.csv")
# head(registros)
# dim(registros)
# #checando NAs----
# names(registros)
# #extrai os valores ambientais em cada registro e inclui o numero do pixel
# spp.env <- extract(predictors,registros[,c("lon","lat")],cellnumbers=T)
# #junta com os registros para ter coordenadas, espécies e ambiente no mesmo
# registros2 <- cbind(registros,spp.env)
#
# head(registros2)
#
# dim(registros2)
#
# #filtra pelo ambiente, os casos (linhas da tabela ambiental que estejam completos)
# spp.env.filt <- registros2[complete.cases(spp.env),]
# head(spp.env.filt)
# dim(spp.env.filt)
#
#
# # agora o passo é tirar as linhas que tem espécie NOVA igual E ambiente igual
# a <- spp.env.filt[,c("sp.final","cells")]
# spp.filt <- spp.env.filt[!duplicated(a),]
#
# head(spp.filt)
# dim(spp.filt)
# ss <- c("c","b","a", "b", "a", "c")
# unique(sort(ss))
# sort(unique(ss))
# table(ss)
# names.sp <- unique(spp.filt$sp.final)
# names.sp <- sort(names.sp)
# spp.filt$sp.final[1:10]
# spp.filt$sp.final[270:280]
#
#
# names.sp <- names.sp[table(spp.filt$sp.final)>10] #considerando apenas spp com mais de 10 registros. AGORA É 10 REGISTROS UNICOS, SEM DUPLICADOS ENTRE PIXEIS
# length(names.sp)###O NUMERO DE ESPECIES CON N>10 MUDOU, CAIU PARA 2115 ANTES ERA 2200 E POUCO --com nomes corrigidos caiu para 2114
#
# occs <- spp.filt[,c("sp.final","lon","lat")]
# occs <- occs[occs$sp.final%in%names.sp,]
# dim(occs)
# names(occs)
#
# rm(spp.filt)
# rm(a)
# rm(registros)
# rm(registros2)
# rm(spp.env)
# rm(spp.env.filt)
# #### ATE AQUI, O RESTO FICOU IGUAL
#
# write.csv(occs,"./data/registroslimpos.csv")
occs <- read.csv("./data/registroslimpos.csv",row.names=1,col.names=c("sp","lon","lat"))
dim(occs)
head(occs)
dim(occs)
is(occs)
class(occs$sp)
names.sp <- unique(occs$sp)
names.sp==sort(names.sp)
names.sp <- sort(names.sp)####ISTO AQUI É ESSENCIAL PARA FAZER O TABLE
names(table(occs$sp))==names.sp
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


