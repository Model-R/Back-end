#16fev2016
#### SEPAREI ESTA PARTE DO SCRIPT POIS NÅO DEVE MAIS SER RODADA MAS É BOM TER.
# INCLUI LIMPEZA TAXONOMICA E FILTROS DE 1. NAS 2. DUPLICADOS 3. REGISTROS UNICOS POR PIXEL 4. ESPECIES COM ESTRITAMENTE MAIS DE 10 OCORRENCIAS

###RODEI ISTO PARA CRIAR A NOVA TABELA DE REGISTROS COM OS NOMES CORRIGIDOS. AGORA PODE PULAR PARA A SEGUINTE PARTE
# RODEI DE NOVO EM 16/02/16 PARA ATUALIZAR A TABELA SEGUNDO A TAXONOMIA - MUDOU MUITO POUCO MAS AGORA SIM FICOU PERFEITO (espero). OS NOMES ESTÅO COM AUTORES E FICARAM COMPLETOS COM SUBESPECIES INCLUSIVE, POIS ISSO ESTAVA DANDO ERRO. são 3918 nomes corretos e 2113 especies com mais de 10 ocorrencias.

#O RESULTADO FOI CRIAR
# 1 registrosnomesfinais.csv
# registroslimpos.csv
#tem também ./data/tax flora/tableoccs.csv com as ocorrencias por especie



predictors <-  raster::stack(list.files("./env",pattern="1K",full.names = T)[1])
####Ler os dados de ocorrencia - tabela por registro----

registros <- read.csv("./data/Atlantic_Domain_Ari_RevTaxon_Final.csv",sep=";")
nomes <- read.delim("./data/tax flora/Lista final.txt")
#
head(registros)
names(registros)

 registros$old.sp <- paste(registros$genero,registros$especie)
 length(registros$old.sp)
 length(unique(registros$old.sp))

dim(nomes)
names(nomes)
length(unique(nomes$FINAL))

library(data.table)
# nomes velhos

setdiff(nomes$original.search,registros$old.sp)
setdiff(registros$old.sp,nomes$original.search)

length(nomes$original.search)
length(unique(nomes$original.search))
length(unique(nomes$FINAL))

length(unique(registros$old.sp))

setdiff(unique(registros$old.sp),unique(nomes$original.search))
setdiff(unique(nomes$original.search),unique(registros$old.sp))
nomes$old.sp <- nomes$original.search
library(dplyr)
registros2 <- left_join(registros,nomes,by="old.sp")
head(registros2)
length(unique(registros2$FINAL))
# write.csv(registros2, file = "./data/1 registrosnomesfinais.csv")
#
# registros <- read.csv("./data/1 registrosnomesfinais.csv")
head(registros)
dim(registros)
# #checando NAs----
names(registros)
# #extrai os valores ambientais em cada registro e inclui o numero do pixel
 spp.env <- extract(predictors,registros[,c("lon","lat")],cellnumbers=T)
# #junta com os registros para ter coordenadas, espécies e ambiente no mesmo
 registros3 <- cbind(registros,spp.env)
#
 head(registros3)
#
 dim(registros3)
#
# #filtra pelo ambiente, os casos (linhas da tabela ambiental que estejam completos)
 spp.env.filt <- registros3[complete.cases(spp.env),]
 head(spp.env.filt)
 dim(spp.env.filt)
#
#
# # agora o passo é tirar as linhas que tem espécie NOVA igual E ambiente igual
 a <- spp.env.filt[,c("FINAL","cells")]
 spp.filt <- spp.env.filt[!duplicated(a),]
 spp.filt$FINAL <- droplevels(spp.filt$FINAL)
 head(spp.filt)
 dim(spp.filt)

 names.sp <- unique(spp.filt$FINAL)

 names.sp <- sort(names.sp)

 names.sp <- names.sp[table(spp.filt$FINAL)>10] #considerando apenas spp com mais de 10 registros. AGORA É 10 REGISTROS UNICOS, SEM DUPLICADOS ENTRE PIXEIS
 length(names.sp)###O NUMERO DE ESPECIES CON N>10 MUDOU, CAIU PARA 2115 ANTES ERA 2200 E POUCO --com nomes corrigidos caiu para 2114


 occs <- spp.filt[,c("FINAL","lon","lat")]
 occs <- occs[occs$FINAL%in%names.sp,]
 dim(occs)
 names(occs)
# write.csv(as.data.frame(table(droplevels(occs$FINAL))),"./data/tax flora/tableoccs.csv")
 View(as.data.frame(table(droplevels(occs$FINAL))))
#write.csv(occs,"./data/registroslimpos.csv")