#Teste da influencia do buffer no TSS
library("raster")
library("XML")
library("snowfall")
library("rJava")

set.seed(712)
N <- sample(1:length(names.sp),30)

args(dismo.mod)

#iniciar snowfall
sfInit(parallel=T,cpus=24)
# 	#exporta variáveis e funções e pacotes
sfExportAll()
sfLibrary(rJava) #removido aspas.
sfLibrary(raster)
sfSource("./fct/modelos.R")
sfSource("./fct/final_model.R")
sfSource("./fct/ensemble.R")

#Com buffer ----
tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = TRUE, buffer.type = "max",
                 seed=712, predictors = predictors, 
                 Bioclim=T, Domain=T, Mahal=T, GLM=T, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "buffermax")

tFinal <- Sys.time()
tFinal - tInicial

#Sem buffer ----
tInicial <- Sys.time()
sfClusterApplyLB(names.sp[N], 
                 dismo.mod, 
                 occs=occs,  
                 buffer = FALSE, buffer.type = "",
                 seed=712, predictors = predictors, 
                 Bioclim=T, Domain=T, Mahal=T, GLM=T, RF=T, SVM=T, SVM2=F, maxent=T, 
                 part=3, n.back=500, output.folder = "sem_buffer")

tFinal <- Sys.time()
tFinal - tInicial

sfStop()
# Analise dos resultados com e sem buffer ----

#organização dos dados

Source("./fct/read.eval1.R")
Source("./fct/read.eval2.R")

#read.eval(sp=names.sp[N][1])
lapply(names.sp[N], read.eval,
       input.folder1="buffermax",
       input.folder2="sem_Buffer",
       output.folder="Evaluate_Flora")

dados<-consolidate.data(stat='mean')

head(dados)
library(ggplot2)
ggplot(dados, aes(x=tratamento, y=tss.buffer.bioclim)) + geom_boxplot(notch=TRUE)
ggplot(dados, aes(x=tratamento, y=tss.buffer.glm)) + geom_boxplot(notch=TRUE)
ggplot(dados, aes(x=tratamento, y=tss.buffer.rf)) + geom_boxplot(notch=TRUE)
ggplot(dados, aes(x=tratamento, y=tss.buffer.svm)) + geom_boxplot(notch=TRUE)
dev.off()

# Sem Snowfall ----

##Com Buffer
lapply(names.sp[N],
       dismo.mod,
       occs = occs,
       predictors = predictors,
       buffer = T,
       buffer.type ="max",
       maxent = T,
       Bioclim = T,
       Domain = T,
       Mahal = T,
       GLM = T,
       RF = T,
       SVM = T,
       SVM2 = F,
       part = 3,
       seed = 712,
       output.folder = "buffermax")

##Sem Buffer
lapply(names.sp[N],
       dismo.mod,
       occs = occs,
       predictors = predictors,
       buffer = F,
       buffer.type ="",
       maxent = T,
       Bioclim = T,
       Domain = T,
       Mahal = T,
       GLM = T,
       RF = T,
       SVM = T,
       SVM2 = F,
       part = 3,
       seed = 712,
       output.folder = "sem_buffer")