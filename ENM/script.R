quit()#leaving R, when using rgedit like me :)
#SSH conection with buriti
ssh felipe@buriti.lncc.br

# creating tmux session
tmux
# IF a tmux session were previously crreated, get attached to it:
tmux attach

# running an empty R session, if necessary
R

# Loading libraries:
library("raster")
library("XML")
library("snowfall")
#install.packages("rjava")
library("rJava")
library("rgdal")

# Loading environmental data, study area mask
predictors <-  stack(list.files("./env",pattern="1K",full.names = T)[1])
# Cortando pela MataAtlantica:
crop <- readOGR(dsn="./data", layer="Bioma_MA1148")
mascara <- rasterize(crop, predictors[[1]], field=crop@data$FID)
mascara <- crop(mascara,crop)

# FLORA ----
# Loading occorrences:
occs <- read.csv("./data/FLORA_occs_final.csv")
head(occs)

# Defining names to be modelled (after taxa and spacial cleaning) see: "./script.limpeza.registros.R"
names.sp <- unique(occs$sp)
length(names.sp)

# set.seed(712)
# N <- sample(1:length(names.sp),length(names.sp))
names.sp[N]
names.sp[N][1:123] # Modelando con n igual a anfibios
names.sp[N][124:151] # Modelando igual Aves
names.sp[N][151:length(names.sp)] # Modelando tudo

# MODELOS.R----
source("./fct/modelos.R")
args(dismo.mod)

#iniciar snowfall
sfInit(parallel=T,cpus=1, slaveOutfile="/home/felipe/FLORA_040316_modelagem.log")
# 	#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")


#Com buffer ----
tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N][1:123], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = TRUE, buffer.type = "max",
                 seed=712, predictors = predictors, 
                 Bioclim=F, Domain=F, Mahal=F, GLM=F, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "FLORA_buffermax_040316", mask=mascara)

tFinal <- Sys.time()
tFinal - tInicial
sfStop()

#Sem buffer ----
#iniciar snowfall
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/FLORA_saida_modelagem_sembuffer.log")

#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")

tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N,], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = FALSE, buffer.type = "",
                 seed=712, predictors = predictors, 
                 Bioclim=T, Domain=T, Mahal=F, GLM=T, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "FLORA_sem_buffer")

tFinal <- Sys.time()
tFinal - tInicial
sfStop()


# Analise do TSS ----
source("./fct/read.eval1.R")

lapply(names.sp, read.eval,
       input.folder1="buffermax",
       input.folder2="sem_buffer",
       output.folder="FLORA_TSS_Evaluate2")

#Criando Grafico
consolidate.data(input.folder1="FLORA_TSS_Evaluate",stat='mean')

# Tabela com n das spp----
species.table("FLORA_TSS_Evaluate",table.name="FLORA.n.sp.modelado")

# Cortando pela MataAtlantica: ----
source('./fct/read.eval1.R')
#args(rasterCrop)
#lapply(names.sp, rasterCrop, input.folder="buffermax", mascara=mascara, crop=crop)
#lapply(names.sp, rasterCrop, input.folder="sem_buffer", mascara=mascara, crop=crop)

# finalModel ----
source("./fct/final_model.R")
args(finalModel)
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/FLORA_040316_finalModel.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/final_model.R")

tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N][1:123],
                 finalModel,
                 input.folder="FLORA_buffermax2")
tFinal <- Sys.time()
tFinal - tInicial
sfStop()


# ENSEMBLE ----
source("./fct/ensemble.R")
args(ensemble)
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/FLORA_040316_ensemble.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/ensemble.R")

tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N][1:123],
                 ensemble,
                 input.folder1="FLORA_buffermax2",occs=occs, consensus=TRUE, consensus.level=0.5)
tFinal <- Sys.time()
tFinal - tInicial
sfStop()

# Analise finalModel() ----
source('./fct/read.eval1.R')
args(partitionEnsemble)
ensemble_list <- list()
ensemble_final <- lapply(names.sp[N,], partitionEnsemble, projeto="FLORA_buffermax2",
                         input.folder='presfinal',
                         algoritmos = c("maxent", "rf", "svm"))

#removendo NA do resultado
ensemble_final <- na.omit(rbindlist(ensemble_final))

library(ggplot2)

#plot diferenciando algoritmos
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  facet_grid(.~algoritmo) +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
#Salvando plot
ggsave(filename="./FLORA_algoritmos.png")

#Gerando plot unico
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
#Salvando plot
ggsave(filename="./FLORA.png")

# Analise ensemble()----
source('./fct/read.eval1.R')
args(spEnsemble)
ensemble_final <- lapply(names.sp[,N], spEnsemble, projeto="FLORA_buffermax2")
#removendo NA do resultado
ensemble_final <- na.omit(rbindlist(ensemble_final))

#gerando plot
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
#salvando plot
ggsave(filename="./FLORA_ensemble.png")

#Salvando raster final----
source("./fct/read.eval1.R")
args(rast.convert)
lapply(names.sp[N,], rast.convert, projeto="FLORA_buffermax2", input.folder="ensemble")

#####
# AVES
#####
# Loading occorrences:
occs <- read.csv("./data/AVES_occs_final.csv",row.names=1,col.names=c("sp","lon","lat"))
# Loading names to be modelled (after taxa and spacial cleaning) see: "./script.limpeza.registros.R"
names.sp <- unique(occs$sp)
length(names.sp)


#selecionando aleatoriamente
set.seed(712)
N <- sample(1:length(names.sp),length(names.sp))
names.sp[N]
names.sp[N][1:123] # Modelando n igual a anfibios
names.sp[N][123:length(names.sp)] # Modelando tudo

#Iniciando SnowFall----
#iniciar snowfall
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/AVES_040316_modelagem.log")

#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")

#AVES Com buffer ----
args(dismo.mod)
tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N][1:123], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = TRUE, buffer.type = "max",
                 seed=712, predictors = predictors, 
                 Bioclim=F, Domain=F, Mahal=F, GLM=F, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "AVES_final", mask=mascara, crop=crop)
tFinal <- Sys.time()
tFinal - tInicial
sfStop()


#Sem buffer ----
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/AVES_saida_semBuffer.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")

tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N,], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = FALSE, buffer.type = "",
                 seed=712, predictors = predictors, 
                 Bioclim=T, Domain=T, Mahal=F, GLM=T, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "AVES_sem_buffer")

tFinal <- Sys.time()
tFinal - tInicial
sfStop()

#Analise do TSS
source("./fct/read.eval1.R")

lapply(names.sp[N,], read.eval,
       input.folder1="AVES_buffermax",
       input.folder2="AVES_sem_buffer",
       output.folder="AVES_TSS_Evaluate")

#Criando Grafico
consolidate.data(input.folder1="AVES_TSS_Evaluate",stat='mean')

# Tabela com n das spp----
species.table("AVES_TSS_Evaluate",table.name="AVES.n.sp.modelado")

# Cortando pela Mata Atlantica
lapply(names.sp[N,][1:30], rasterCrop, input.folder="AVES_buffermax", mascara=mascara, crop=crop)
lapply(names.sp[N,][1:30], rasterCrop, input.folder="AVES_sem_buffer", mascara=mascara, crop=crop)

# finalModel ----
args(finalModel)
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/AVES_040316_finalModel.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/final_model.R")

tInicial <- Sys.time()
lapply(sort(names.sp[N][1:123]),
       finalModel,
       input.folder="AVES_final")
tFinal <- Sys.time()
tFinal - tInicial

# ENSEMBLE ----
args(ensemble)
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/ANFIBIOS_040316_ensemble.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/ensemble.R")

tInicial <- Sys.time()
sfClusterApplyLB(sort(names.sp[N,][1:123]),
                 ensemble,
                 input.folder1="AVES_final",occs=occs, consensus=TRUE, consensus.level=0.5)

tFinal <- Sys.time()
tFinal - tInicial
sfStop()

# Analise finalModel ----
source('./fct/read.eval1.R')
ensemble_list <- list()
ensemble_final <- lapply(names.sp[N,], partitionEnsemble, projeto="./AVES_buffermax2",
                         input.folder='presfinal',
                         algoritmos = c("maxent", "rf", "svm")
)

ensemble_final <- na.omit(rbindlist(ensemble_final))

library(ggplot2)
#plot diferenciando algoritmos
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  facet_grid(.~algoritmo) +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
ggsave(filename="./AVES_algoritmos.png")

#Salvando plot unico
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
ggsave(filename="./AVES.png")

# Analise ensemble() ----
source('./fct/read.eval1.R')
args(spEnsemble)
ensemble_final <- lapply(names.sp[N,], spEnsemble, projeto="AVES_buffermax2")
#removendo NA do resultado
ensemble_final <- na.omit(rbindlist(ensemble_final))

#gerando plot
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
#salvando plot
ggsave(filename="./AVES_ensemble.png")

#Salvando raster final ----
source("./fct/read.eval1.R")
args(rast.convert)
lapply(names.sp[N,], rast.convert, projeto="AVES_buffermax2", input.folder="ensemble")

#####
#ANFIBIOS
#####

# Loading occorrences:
occs <- read.csv("./data/ANFIBIOS_occs_final.csv",row.names=1,col.names=c("sp","lon","lat"))
# Loading names to be modelled (after taxa and spacial cleaning) see: "./script.limpeza.registros.R"
names.sp <- unique(occs$sp)
length(names.sp)

#selecionando aleatoriamente
set.seed(712)
N <- sample(1:length(names.sp),length(names.sp))
names.sp[N]

#Iniciando SnowFall----
#iniciar snowfall
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/ANFIBIOS_040316_modelagem.log")

#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")

#ANFIBIOS Com buffer ----
tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = TRUE, buffer.type = "max",
                 seed=712, predictors = predictors, 
                 Bioclim=F, Domain=F, Mahal=F, GLM=F, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "ANFIBIOS_final", mask=mascara, crop=crop)
tFinal <- Sys.time()
tFinal - tInicial
sfStop()

#Sem  Buffer-----
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/ANFIBIOS_saida_sembuffer_modelagem.log")

#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")

#ANFIBIOS sem buffer 
tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N,], 
                 dismo.mod,
                 occs=occs,  
                 buffer = FALSE, buffer.type = "",
                 seed=712, predictors = predictors, 
                 Bioclim=T, Domain=T, Mahal=F, GLM=T, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "ANFIBIOS_sem_buffer")

tFinal <- Sys.time()
tFinal - tInicial
sfStop()

#Analise do TSS ----
source("./fct/read.eval1.R")

lapply(names.sp[N,], read.eval,
       input.folder1="ANFIBIOS_buffermax",
       input.folder2="ANFIBIOS_sem_buffer",
       output.folder="ANFIBIOS_TSS_Evaluate")

#Criando Grafico
consolidate.data(input.folder1="ANFIBIOS_TSS_Evaluate",stat='mean')

# Tabela com n das spp----
species.table("ANFIBIOS_TSS_Evaluate",table.name="ANFIBIOS.n.sp.modelado")

# Cortando pela Mata Atlantica
sfInit(parallel=T,cpus=5, slaveOutfile="/home/felipe/ANFIBIOS_saida_rasterCrop.log")
sfClusterApplyLB(names.sp[N,], rasterCrop, input.folder="ANFIBIOS_buffermax", mascara=mascara, crop=crop)
sfStop()
lapply(names.sp[N][1:30], rasterCrop, input.folder="ANFIBIOS_sem_buffer", mascara=mascara, crop=crop)

# finalModel ----
source('./fct/final_model.R')
args(finalModel)
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/ANFIBIOS_040316_finalModel.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/final_model.R")

tInicial <- Sys.time()
sfClusterApplyLB(sort(names.sp[N]),
                 finalModel,
                 input.folder="ANFIBIOS_final")

tFinal <- Sys.time()
tFinal - tInicial
sfStop()

# ENSEMBLE ----
args(ensemble)
sfInit(parallel=T,cpus=10, slaveOutfile="/home/felipe/ANFIBIOS_040316_ensemble.log")
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/ensemble.R")

tInicial <- Sys.time()
sfClusterApplyLB(sort(names.sp[N]),
                 ensemble,
                 input.folder1="ANFIBIOS_final",occs=occs, consensus=TRUE, consensus.level=0.5)

tFinal <- Sys.time()
tFinal - tInicial
sfStop()

# Analise finalModel ----
source('./fct/read.eval1.R')
ensemble_list <- list()
ensemble_final <- lapply(names.sp[N,], partitionEnsemble, projeto="./ANFIBIOS_buffermax2",
                         input.folder='presfinal',
                         algoritmos = c("maxent", "rf", "svm")
)

ensemble_final <- na.omit(rbindlist(ensemble_final))

library(ggplot2)
#plot diferenciando algoritmos
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  facet_grid(.~algoritmo) +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
ggsave(filename="./ANFIBIOS_algoritmos.png")

#Salvando plot unico
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
ggsave(filename="./ANFIBIOS.png")

# Analise ensemble() -----
source('./fct/read.eval1.R')
args(spEnsemble)
ensemble_final <- lapply(names.sp[N,], spEnsemble, projeto="ANFIBIOS_buffermax2")
#removendo NA do resultado
ensemble_final <- na.omit(rbindlist(ensemble_final))

#gerando plot
ggplot(ensemble_final, aes(x=n.sp, y=ensemble, colour=ensemble)) + geom_point() +  scale_colour_gradientn(colours = c("orange", "red", "darkred")) + labs(x = "n", y = "% Consenso entre partições") + theme(text=element_text(size=16))
#salvando plot
ggsave(filename="./ANFIBIOS_ensemble.png")

#Salvando raster final ----
source("./fct/read.eval1.R")
args(rast.convert)
lapply(names.sp[N,], rast.convert, projeto="ANFIBIOS_buffermax2", input.folder="ensemble")

#######################
#vendo quantidade de arquivos na pasta:
for especie in *; do echo "$especie" $(ls -1 "$especie" | wc -l); done
for especie in *; do echo "$especie" $(ls -1 "$especie" | wc -l); done | grep -v 121

#########################
#compactando os dados
tar cjf NomeArquivoCompactado.tar.bz2 NomeDaPastaASerCompactada

#DESCOMPACTANDO
tar -xvfj NomeArquivoCompactado.tar.bz2

#fazendo download 
#Fora do buriti executar
sftp felipe@buriti.lncc.br
get nome_do_arquivo
#ou
get -r nome_da_pasta

######
#Problema com RJava: Rodar no bash
R CMD javareconf
