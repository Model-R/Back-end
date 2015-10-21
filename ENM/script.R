#carrega a função dismo.mod ----
source("./fct/modelos.R")
library("raster")
library("XML")
library("snowfall")
library("rJava")
#Ler os raster ambientais----
file <-  list.files("./env",full.names = T)[1]
predictors <- stack(file)
plot(predictors)

####Ler os dados de ocorrencia - tabela por registro----
registros<-read.csv("./data/Atlantic_Domain_Ari_RevTaxon_Final.csv",sep=";")
head(registros)
#checando NAs----
ex <- extract(predictors,registros[,6:5])
spp.filt <- registros[complete.cases(ex),]

dim(registros)
dim(spp.filt)
names(spp.filt)
spp.filt$sp <- paste(spp.filt$genero,spp.filt$especie)
occs <- spp.filt
args(dismo.mod)
names.sp <- unique(spp.filt$sp)
names <- names.sp[table(spp.filt$sp)>10]
dismo.mod(names[3],predictors = predictors,Bioclim = T,maxent = T,occs=spp.filt,seed=123)
}sp <- names.sp[1]
#iniciar snowfall 
sfInit(parallel=T,cpus=24)
	#exporta variáveis e funções e pacotes
	sfExportAll()
	sfLibrary("rJava")
	sfSource("./fct/modelos.R")
	
	#executa a função 
	sfClusterApplyLB(names.sp,dismo.mod)

sfStop()

q()