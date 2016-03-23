# Script for Flora's taxonomic cleaning
# and Flora, anphibian and birds data organization:
# (only species with more than 10 unique occurrences and atleast 50% of occurrences in Altantic Rain Forest will be modeled) -> spactial cleaning

# As Results:
# registrosnomesfinais.csv # Taxonomic cleanning
# tableoccs.csv com as ocorrencias por especie #
# registroslimpos.csv # Spacial data cleaning
# Flora_occs_final.csv # Spacial data cleaning II

# Library required
library(raster)
library(rgdal)
library(data.table)
library(dplyr)

# environmental and spatial data
predictors <-  raster::stack(list.files("./env",pattern="1K",full.names = T)[1])

# Altantic Rain forest limit:
crop <- readOGR(dsn="./data", layer="Bioma_MA1148")
mascara<-rasterize(crop, predictors[[1]], field=crop@data$FID)
#zonal_maks <- crop(mascara,crop)

# FLORA:
## Taxonomic cleaning:----
registros <- read.csv("./data/Atlantic_Domain_Ari_RevTaxon_Final.csv",sep=";")
nomes <- read.delim("./data/tax flora/lista final.txt")

### check results
head(registros)
names(registros)

### combining genero and species name
registros$old.sp <- paste(registros$genero,registros$especie)
length(registros$old.sp)
length(unique(registros$old.sp))

dim(nomes)
names(nomes)
length(unique(nomes$FINAL))

### nomes velhos
#setdiff(nomes$original.search,registros$old.sp)
#setdiff(registros$old.sp,nomes$original.search)
#length(nomes$original.search)
#length(unique(nomes$original.search))
#length(unique(nomes$FINAL))
#length(unique(registros$old.sp))
#setdiff(unique(registros$old.sp),unique(nomes$original.search))
#setdiff(unique(nomes$original.search),unique(registros$old.sp))
nomes$old.sp <- nomes$original.search

registros2 <- left_join(registros,nomes,by="old.sp")
head(registros2)
length(unique(registros2$FINAL))

### Taxonomic cleaning concluded: saving result:
# write.csv(registros2, file = "./data/1 registrosnomesfinais.csv")
## Spatial data cleaning:----
# registros <- read.csv("./data/1 registrosnomesfinais.csv")
head(registros)
dim(registros)

### checking NAs (occurrences outside the study area)
names(registros)

### extracting values from predictor raster layer and pixel number
spp.env <- extract(predictors,registros[,c("lon","lat")],cellnumbers=T)

### Combinig extract with occurrence data
registros2 <- cbind(registros,spp.env)

head(registros2)
dim(registros2)

### filtering data with complete environmental values
spp.env.filt <- registros2[complete.cases(spp.env),]

head(spp.env.filt)
dim(spp.env.filt)

### Removing occurrences with duplicated points (tow points of spp in same pixel)
a <- spp.env.filt[,c("FINAL","cells")]
spp.filt <- spp.env.filt[!duplicated(a),]
spp.filt$FINAL <- droplevels(spp.filt$FINAL)

head(spp.filt)
dim(spp.filt)

### Saving species names
names.sp <- unique(spp.filt$FINAL)
names.sp <- sort(names.sp)

### Filtering species with more than 10 unique occurrences
names.sp <- names.sp[table(spp.filt$FINAL)>10]

length(names.sp)###O NUMERO DE ESPECIES CON N>10 MUDOU, CAIU PARA 2115 ANTES ERA 2200 E POUCO --com nomes corrigidos caiu para 2114

### Organizing occurrences (occs) data with species that completed all requirements
occs <- spp.filt[,c("FINAL","lon","lat")]
occs <- occs[occs$FINAL%in%names.sp,]
dim(occs)
names(occs)
### Spatial data cleaning concluded: saving results:
# write.csv(as.data.frame(table(droplevels(occs$FINAL))),"./data/tax flora/tableoccs.csv")
#View(as.data.frame(table(droplevels(occs$FINAL))))
#write.csv(occs,"./data/registroslimpos.csv")
## Spatial data cleaning II: removing species with marginal distribution on Altantic Rainforest ----
#occs <- read.csv("./data/registroslimpos.csv",row.names=1,col.names=c("sp","lon","lat"))
dim(occs)
names.sp <- unique(occs$sp)
names.sp <- sort(names.sp)####ISTO AQUI É ESSENCIAL PARA FAZER O TABLE
length(names.sp)
 
### Filtering species with marginal distribution on Mata Atlantica
table <- as.data.frame(table(occs$sp)) #table with occurences by spp

# Extracting values from Atlantic Rainforest mask to identify which occurrences are inside of it
spp.env2 <- extract(mascara,occs[,c("lon","lat")],cellnumbers=T)

spp.env3 <- cbind(occs,spp.env2)
head(spp.env3)
spp.env3 <- table(spp.env3$sp, spp.env3$layer)
spp.env3 <- as.data.frame(spp.env3)
spp.env3 <- spp.env3[,c(1,3)]

### Making a table with freq per species of occs in Atlantic Rain Forest
table2 <- merge(table, spp.env3, by.x="Var1", by.y="Var1")
head(table2)
colnames(table2) <- c('sp','n','n.ma')

### Estimating a percentual of occs in Atlantic Rainforest from all spp occs
table2$per <- (table2$n.ma/table2$n)
head(table2)

### Keeping onle species with more than 50% of occs in Atlantic Rainforest 
table3 <- table2[which(table2$per>=.5),]
head(table3)

### Updating names.sp with spp names that comply with all procedures 
names.sp <- unique(table3$sp)
length(names.sp)

### Updating occs table with occs of spp that comply all procedures 
occs_final <- occs[occs$sp %in% names.sp,]
dim(occs_final)
unique(occs_final$sp)
# Spatial data cleaning II done. Saving results
#write.csv(occs_final,"./data/FLORA_occs_final.csv")
 
# ANFIBIOS: ----
registros <- read.csv("./data/anfibios.csv")
head(registros)
registros$sp <- paste(registros$Especies)
 
length(unique(registros$sp))
 
### Spatial datacleaning: checking NAs ----
#### extracting environmental values and pixel number for each occs
spp.env <- extract(predictors,registros[,c("lon","lat")],cellnumbers=T)

### Combinig extract with occurrence data
registros2 <- cbind(registros,spp.env)
head(registros2)
dim(registros2)
 
### filtering data with complete environmental values
spp.env.filt <- registros2[complete.cases(spp.env),]
head(spp.env.filt)
dim(spp.env.filt)
 
### Removing occurrences with duplicated points (tow points of spp in same pixel)
a <- spp.env.filt[,c("sp","cells")]
spp.filt <- spp.env.filt[!duplicated(a),]
dim(spp.filt)

### Saving species names
names.sp <- unique(spp.filt$sp)
names.sp <- sort(names.sp)
head(names.sp)
length(names.sp)

### Filtering species with more than 10 unique occurrences
names.sp <- names.sp[table(spp.filt$sp)>10]
length(names.sp)

### Organizing occurrences (occs) data with species that completed all requirements
occs <- spp.filt[,c("sp","lon","lat")]
occs <- occs[occs$sp %in% names.sp,]
dim(occs)
 
rm(spp.filt)
rm(a)
rm(registros)
rm(registros2)
rm(spp.env)
rm(spp.env.filt)
 
### Spatial data cleaning concluded: saving results:
# write.csv(occs,"./data/ANFIBIOS_registroslimpos.csv")
### Spatial data cleaning II: removing species withmarginal distribution on Altantic Rain Forest ----
# occs <- read.csv("./data/ANFIBIOS_registroslimpos.csv",row.names=1,col.names=c("sp","lon","lat"))
names.sp <- unique(occs$sp)
names.sp <- sort(names.sp)####ISTO AQUI É ESSENCIAL PARA FAZER O TABLE
length(names.sp)

### Filtering species with marginal distribution on Mata Atlantica
table <- as.data.frame(table(occs$sp))
head(table)
tail(table)

# Extracting values from Atlantic Rainforest mask to identify which occurrences are inside of it
spp.env2 <- extract(mascara,occs[,c("lon","lat")],cellnumbers=T)
head(spp.env2)

spp.env3 <- cbind(occs,spp.env2)
spp.env3 <- table(spp.env3$sp, spp.env3$layer)
head(spp.env3)
colnames(spp.env3)
spp.env3<-as.data.frame(spp.env3)
spp.env3<-spp.env3[,c(1,3)]
head(spp.env3)

### Making a tablewith freq per species of occs in Atlantic Rain Forest
table2 <- merge(table, spp.env3, by.x="Var1", by.y="Var1")
head(table2)
colnames(table2) <- c('sp','n','n.ma')

### Estimating a percentual of occs in Atlantic Rainforest from all spp occs
table2$per <- (table2$n.ma/table2$n)
head(table2)

### Keeping onle species with more than 50% of occs in Atlantic Rainforest 
table3 <- table2[which(table2$per>=.5),]
head(table3)
which(table3$per<.5)	

### Updating names.sp with spp names that comply with all procedures 
names.sp <- unique(table3$sp)
length(names.sp)

### Updating occs table with occs of spp that comply all procedures 
occs_final <- occs[occs$sp %in% names.sp,]
dim(occs_final)

# Spatial data cleaning II done. Saving results
# write.csv(occs_final,"./data/ANFIBIOS_occs_final.csv")
# AVES: ----
occs <- read.csv("./data/AVES_registroslimpos.csv",row.names=1,col.names=c("sp","lon","lat"))

dim(occs)
head(occs)
dim(occs)
### Spatial datacleaning: checking NAs ----
#### extracting environmental values and pixel number for each occs
spp.env <- extract(predictors, occs[,c("lon","lat")], cellnumbers=T)

### Combinig extract with occurrence data
registros2 <- cbind(occs, spp.env)
head(registros2)
dim(registros2)

### filtering data with complete environmental values
spp.env.filt <- registros2[complete.cases(spp.env),]
head(spp.env.filt)
dim(spp.env.filt)

### Removing occurrences with duplicated points (tow points of spp in same pixel)
a <- spp.env.filt[,c("sp","cells")]
spp.filt <- spp.env.filt[!duplicated(a),]
dim(spp.filt)

### Saving species names
names.sp <- unique(spp.filt$sp)
names.sp <- sort(names.sp)
head(names.sp)
length(names.sp)

### Filtering species with more than 10 unique occurrences
names.sp <- droplevels(names.sp[table(spp.filt$sp)>10])
length(names.sp)

### Organizing occurrences (occs) data with species that completed all requirements
occs <- spp.filt[,c("sp","lon","lat")]

occs <- spp.filt[(spp.filt$sp %in% names.sp),c("sp","lon","lat")]
dim(occs)

### Spatial data cleaning concluded: saving results:
# write.csv(occs,"./data/AVES_registroslimpos.csv")
### Spatial data cleaning II: removing species withmarginal distribution on Altantic Rain Forest ----
# occs <- read.csv("./data/AVES_registroslimpos.csv",row.names=1,col.names=c("sp","lon","lat"))
# names.sp <- unique(occs$sp)
# names.sp <- sort(names.sp)####ISTO AQUI É ESSENCIAL PARA FAZER O TABLE
length(names.sp)

### Filtering species with marginal distribution on Mata Atlantica
table <- as.data.frame(table(occs$sp))
head(table)
tail(table)

# Extracting values from Atlantic Rainforest mask to identify which occurrences are inside of it
spp.env2 <- extract(mascara,occs[,c("lon","lat")],cellnumbers=T)
head(spp.env2)

spp.env3 <- cbind(occs,spp.env2)
spp.env3 <- table(spp.env3$sp, spp.env3$layer)
head(spp.env3)
colnames(spp.env3)
spp.env3<-as.data.frame(spp.env3)
spp.env3<-spp.env3[,c(1,3)]
head(spp.env3)

### Making a tablewith freq per species of occs in Atlantic Rain Forest
table2 <- merge(table, spp.env3, by.x="Var1", by.y="Var1")
head(table2)
colnames(table2) <- c('sp','n','n.ma')

### Estimating a percentual of occs in Atlantic Rainforest from all spp occs
table2$per <- (table2$n.ma/table2$n)
head(table2)
tail(table2)

### Keeping onle species with more than 50% of occs in Atlantic Rainforest 
table3 <- table2[which(table2$per>=.5),]
head(table3)
which(table3$per<.5)	

### Updating names.sp with spp names that comply with all procedures 
names.sp <- droplevels(unique(table3$sp))
length(names.sp)

### Updating occs table with occs of spp that comply all procedures 
occs_final <- occs[occs$sp %in% names.sp,]
dim(occs_final)

# Spatial data cleaning II done. Saving results
# write.csv(occs_final,"./data/AVES_occs_final.csv")